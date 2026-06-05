#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.samplesheet = "samplesheet.csv"
params.outdir = null
params.kmer_fasta = "${projectDir}/data/all_tagging_kmers.fasta.gz"
params.model_directory = "${projectDir}/models/"
params.help = false

// Help message
def helpMessage() {
    log.info"""
    Usage:
        nextflow run horoscope.nf --samplesheet <path> --outdir <path>

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
      - NORMALIZATION:        Normalization mode ('global', 'chromosome', 'p_arm', 'q_arm', or a float value)
    """.stripIndent()
}

// ---------------------------------------------------------------------------
// Processes
// ---------------------------------------------------------------------------

process DBG_LIGHTER {
    tag "LIGHTER_${name}"

    input:
    tuple val(name), path(input_file), val(file_type)

    output:
    tuple val(name), path("${name}.cor.fastq.gz")

    script:
    """
    if [[ "${file_type}" != "fastq.gz" && "${file_type}" != "fastq" ]]; then
        echo "Read correction with lighter only supports FILE_TYPE='fastq.gz' or 'fastq' for sample ${name}; got '${file_type}'" >&2
        exit 1
    fi

    if [[ "${input_file}" != *.fastq.gz && "${input_file}" != *.fq.gz ]]; then
        echo "Read correction with lighter requires a .fastq.gz or .fq.gz input file for sample ${name}; got '${input_file}'" >&2
        exit 1
    fi

    lighter -t ${task.cpus} -r "${input_file}" -trim -discard -k 23 3100000000 0.188
    """
}

process DBG_BCALM {
    tag "BCALM_${name}"

    input:
    tuple val(name), path(corrected_fastq)

    output:
    tuple val(name), path("${name}.unitigs.fa.gz")

    script:
    """
    echo "dbg_bcalm placeholder: not implemented yet for sample ${name}" >&2
    exit 1
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

    input:
    tuple val(name), path(input_file)
    path(kmer_fasta)
    path(dbg_script)

    output:
    tuple val(name), path("${name}.kmer.tsv")

    script:
    """
    ### rewrite fasta to a list of kmers
    cat "${kmer_fasta}" | awk 'NR % 2 == 0' > kmer_list.txt

    ### create local uncompressed copy of the dbg file in task work directory
    gzip -dc "${input_file}" > dbg.uncompressed.fa

    ### subset dbg to those contigs with one matching kmers and write to new file
    pv dbg.uncompressed.fa | grep -B 1 -f kmer_list.txt > "${name}".kmers.zgrep.txt
    

    ### call python file to convert the zgrep output into a kmer table
    python3 ${dbg_script} \\
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

    input:
    tuple val(name), path(input_file), val(file_type), val(cram_ref), path(kmer_fasta)

    output:
    tuple val(name), path("${name}.kmer.tsv")

    script:
    """
    if [[ "${file_type}" == "fastq" || "${file_type}" == "fasta" ]]; then
        if [[ "${input_file}" == *.gz || "${input_file}" == *.bgz ]]; then
            jellyfish count -m 61 -s 2G -t ${task.cpus} --if "${kmer_fasta}" -C <(zcat "${input_file}") -o "${name}.jf"
        else
            jellyfish count -m 61 -s 2G -t ${task.cpus} --if "${kmer_fasta}" -C <(cat "${input_file}") -o "${name}.jf"
        fi
    elif [[ "${file_type}" == "bam" ]]; then
        jellyfish count -m 61 -s 2G -t ${task.cpus} --if "${kmer_fasta}" -C <(samtools fasta "${input_file}") -o "${name}.jf"
    elif [[ "${file_type}" == "cram" ]]; then
        if [[ "${cram_ref}" == "NA" || -z "${cram_ref}" ]]; then
            echo "CRAM input requires CRAM_REFERENCE_PATH for sample ${name}" >&2
            exit 1
        fi
        jellyfish count -m 61 -s 2G -t ${task.cpus} --if "${kmer_fasta}" -C <(samtools fasta --reference "${cram_ref}" "${input_file}") -o "${name}.jf"
    else
        echo "Unsupported FILE_TYPE for GATHER_KMERS_JF: ${file_type}" >&2
        exit 1
    fi

    jellyfish query \\
        "${name}.jf" \\
        -s "${kmer_fasta}" \\
        -o "${name}.kmer.tsv"

    rm "${name}.jf"
    """
}

process FORMAT_AND_MERGE_KMERS {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(kmer_files)
    path(reformat_script)
    path(kmer_info_file)

    output:
    path("final_kmer_merged.tsv")

    script:
    """
    ls -1 ${kmer_files} > all_kmer_files.txt
    python3 ${reformat_script} \\
        --kmer_file_list all_kmer_files.txt \\
        --kmer_info_file ${kmer_info_file}
    """
}

process GENOTYPE {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(merged_kmers)
    path(samplesheet)
    path(genotype_script)

    output:
    path("normalization_metrics.tsv")
    path("centromere_genotyping.tsv")

    script:
    """
    python ${genotype_script} \\
        --kmer_table ${merged_kmers} \\
        --samplesheet ${samplesheet} \\
        --model_directory ${params.model_directory}
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
            def normalizationRaw = (row.NORMALOIZATION ?: row.NORMALIZATION ?: '').toString().trim()
            def normalizationLower = normalizationRaw.toLowerCase()
            def isFloatNormalization = normalizationRaw.isNumber()
            if (!(normalizationLower in ['global', 'chromosome', 'p_arm', 'q_arm']) && !isFloatNormalization) {
                error "Invalid NORMALIZATION value '${normalizationRaw}' for sample ${row.NAME}. Allowed values are 'global', 'chromosome', 'p_arm', 'q_arm', or a float."
            }

            tuple(
                row.NAME,
                file(row.FILE_PATH),
                row.FILE_TYPE,
                row.CRAM_REFERENCE_PATH,
                row.MAKE_DBG.toLowerCase() == 'true',
                normalizationRaw
            )
        }
        .branch {
            _name, _input_file, file_type, _cram_ref, make_dbg, _normalization ->
            to_build_dbg:  make_dbg == true
            existing_dbg:  make_dbg == false && file_type == 'dbg'
            other:         make_dbg == false && file_type != 'dbg'
        }
        .set { branched }

    // Unpack tagging kmers once and reuse in all jellyfish tasks
    UNPACK_KMER_FASTA(channel.value(file(params.kmer_fasta, checkIfExists: true)))

    // Build de Bruijn graphs where requested
    DBG_LIGHTER(
        branched.to_build_dbg.map { row ->
            tuple(row[0], row[1], row[2])
        }
    )

    DBG_BCALM(
        DBG_LIGHTER.out
    )

    // Rows that already have a pre-built DBG file go straight to GATHER_KMERS_DBG
    existing_dbg_for_gather_ch = branched.existing_dbg
        .map { name, input_file, _file_type, _cram_ref, _make_dbg, _normalization ->
            tuple(name, input_file)
        }

    // Combine freshly built DBGs and pre-existing DBG files, then gather k-mers
    GATHER_KMERS_DBG(
        DBG_BCALM.out.mix(existing_dbg_for_gather_ch),
        UNPACK_KMER_FASTA.out,
        file("${projectDir}/scripts/dbg_to_kmer_table.py")
    )

    // All remaining file types go to the Jellyfish-based k-mer counter
    GATHER_KMERS_JF(
        branched.other
            .map { name, input_file, file_type, cram_ref, _make_dbg, _normalization ->
                tuple(name, input_file, file_type, cram_ref)
            }
            .combine(UNPACK_KMER_FASTA.out)
            .map { name, input_file, file_type, cram_ref, unpacked_kmer_fasta ->
                tuple(name, input_file, file_type, cram_ref, unpacked_kmer_fasta)
            }
    )

    // Collect names of all kmer output files into a single manifest
    FORMAT_AND_MERGE_KMERS(
        GATHER_KMERS_DBG.out.mix(GATHER_KMERS_JF.out)
            .map { _name, kmer_file -> kmer_file }
            .collect(),
        file("${projectDir}/scripts/reformat_merge_kmers.py"),
        file("${projectDir}/data/all_tagging_kmers.tsv.gz")
    )

    GENOTYPE(
        FORMAT_AND_MERGE_KMERS.out,
        channel.value(file(params.samplesheet, checkIfExists: true)),
        file("${projectDir}/scripts/genotype.py")
    )
}
