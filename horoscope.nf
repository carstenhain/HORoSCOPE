#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.samplesheet = "samplesheet.csv"
params.outdir = null
params.kmer_fasta = "${projectDir}/data/all_tagging_kmers.fasta.gz"

// Help message
def helpMessage() {
    log.info"""
    Usage:
        nextflow run gather_kmers.nf --samplesheet <path> --outdir <path>

    Required arguments:
        --samplesheet PATH    Path to the semicolon-separated samplesheet file
        --outdir PATH         Directory where result files are stored

    Optional arguments:
        --help                Show this help message and exit

    The samplesheet should be a semicolon-separated file with columns:
      - NAME:                 Sample identifier (string)
      - FILE_PATH:            Path to the input file (string)
      - FILE_TYPE:            Type of the input file (string, e.g. 'bam', 'cram', 'fastq', 'dbg')
      - CRAM_REFERENCE_PATH:  Path to reference FASTA for CRAM decoding (or NA)
      - MAKE_DBG:             Whether to build a de Bruijn graph for this sample (true/false)
    """.stripIndent()
}

// ---------------------------------------------------------------------------
// Processes
// ---------------------------------------------------------------------------

process MAKE_DBG {
    tag "${name}"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(name), path(input_file), val(file_type), val(cram_ref)

    output:
    tuple val(name), path("${name}.unitigs.fa.gz")

    script:
    """
    touch ${name}.unitigs.fa.gz
    """
}

process UNPACK_KMER_FASTA {
    input:
    path(kmer_fasta_gz)

    output:
    path("all_tagging_kmers.fasta")

    script:
    """
    gzip -dc "${kmer_fasta_gz}" > all_tagging_kmers.fasta
    """
}

process GATHER_KMERS_DBG {
    tag "${name}"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(name), path(input_file)

    output:
    tuple val(name), path("${name}.kmer.tsv")

    script:
    """
    ### rewrite fasta to a list of kmers
    zcat ${projectDir}/data/all_tagging_kmers.fasta.gz | awk 'NR % 2 == 0' > kmer_list.txt

    ### create local uncompressed copy of the dbg file in task work directory
    gzip -dc "${input_file}" > dbg.uncompressed.fa

    ### subset dbg to those contigs with one matching kmers and write to new file
    pv dbg.uncompressed.fa | grep -B 1 -f kmer_list.txt > "${name}".kmers.zgrep.txt
    

    ### call python file to convert the zgrep output into a kmer table
    python3 ${projectDir}/scripts/dbg_to_kmer_table.py \\
        --name "${name}" \\
        --zgrep_file "${name}.kmers.zgrep.txt"  \\
        --kmer_list kmer_list.txt \\
        --output "${name}.kmer.tsv"

    rm kmer_list.txt
    rm dbg.uncompressed.fa
    """
}

process GATHER_KMERS_JF {
    tag "${name}"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(name), path(input_file), val(file_type), val(cram_ref), path(kmer_fasta)

    output:
    tuple val(name), path("${name}.kmer.tsv")

    script:
    """
    if [[ "${file_type}" == "fastq" || "${file_type}" == "fasta" ]]; then
        if [[ "${input_file}" == *.gz || "${input_file}" == *.bgz ]]; then
            gzip -dc "${input_file}" | jellyfish count -m 61 -s 2G --if ${kmer_fasta} -t ${task.cpus} -C /dev/stdin -o "${name}.jf"
        else
            cat "${input_file}" | jellyfish count -m 61 -s 2G --if ${kmer_fasta} -t ${task.cpus} -C /dev/stdin -o "${name}.jf"
        fi
    elif [[ "${file_type}" == "bam" ]]; then
        samtools fasta "${input_file}" | jellyfish count -m 61 -s 2G --if ${kmer_fasta} -t ${task.cpus} -C /dev/stdin -o "${name}.jf"
    elif [[ "${file_type}" == "cram" ]]; then
        if [[ "${cram_ref}" == "NA" || -z "${cram_ref}" ]]; then
            echo "CRAM input requires CRAM_REFERENCE_PATH for sample ${name}" >&2
            exit 1
        fi
        samtools fasta -T "${cram_ref}" "${input_file}" | jellyfish count -m 61 -s 2G --if ${kmer_fasta} -t ${task.cpus} -C /dev/stdin -o "${name}.jf"
    else
        echo "Unsupported FILE_TYPE for GATHER_KMERS_JF: ${file_type}" >&2
        exit 1
    fi

    jellyfish query \\
        "${name}.jf" \\
        -s "${kmer_fasta}" \\
        -o "${name}.kmer.tsv"
    """
}

process COLLECT_KMER_OUTPUTS {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(kmer_files)

    output:
    path("all_kmer_files.txt")

    script:
    """
    ls -1 ${kmer_files} > all_kmer_files.txt
    """
}

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------

workflow {

    if (params.help) {
        helpMessage()
        return
    }

    // Parse the semicolon-delimited samplesheet, retaining the MAKE_DBG flag for routing
    channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true, sep: ';')
        .map { row ->
            tuple(
                row.NAME,
                file(row.FILE_PATH),
                row.FILE_TYPE,
                row.CRAM_REFERENCE_PATH,
                row.MAKE_DBG.toLowerCase() == 'true'
            )
        }
        .branch {
            name, input_file, file_type, cram_ref, make_dbg ->
            to_build_dbg:  make_dbg == true
            existing_dbg:  make_dbg == false && file_type == 'dbg'
            other:         make_dbg == false && file_type != 'dbg'
        }
        .set { branched }

    // Unpack tagging kmers once and reuse in all jellyfish tasks
    UNPACK_KMER_FASTA(channel.value(file(params.kmer_fasta, checkIfExists: true)))

    // Build de Bruijn graphs where requested
    MAKE_DBG(
        branched.to_build_dbg.map { name, input_file, file_type, cram_ref, make_dbg ->
            tuple(name, input_file, file_type, cram_ref)
        }
    )

    // Rows that already have a pre-built DBG file go straight to GATHER_KMERS_DBG
    existing_dbg_for_gather_ch = branched.existing_dbg
        .map { name, input_file, file_type, cram_ref, make_dbg ->
            tuple(name, input_file)
        }

    // Combine freshly built DBGs and pre-existing DBG files, then gather k-mers
    GATHER_KMERS_DBG(
        MAKE_DBG.out.mix(existing_dbg_for_gather_ch)
    )

    // All remaining file types go to the Jellyfish-based k-mer counter
    GATHER_KMERS_JF(
        branched.other
            .map { name, input_file, file_type, cram_ref, _make_dbg ->
                tuple(name, input_file, file_type, cram_ref)
            }
            .combine(UNPACK_KMER_FASTA.out)
            .map { name, input_file, file_type, cram_ref, unpacked_kmer_fasta ->
                tuple(name, input_file, file_type, cram_ref, unpacked_kmer_fasta)
            }
    )

    // Collect names of all kmer output files into a single manifest
    COLLECT_KMER_OUTPUTS(
        GATHER_KMERS_DBG.out.mix(GATHER_KMERS_JF.out)
            .map { name, kmer_file -> kmer_file }
            .collect()
    )
}
