# centromere_genotyping
Workflow to genotype the human centromere (architecture and mean alpha-satellite length) from short-read data

ToDo
Command
nextflow run /g/korbel/hain/centromere_genotyping_tool/centromere_genotyping/genotype_centromere.nf -c /g/korbel/hain/centromere_genotyping_tool/centromere_genotyping/nextflow.config --samplesheet /g/korbel/hain/centromere_genotyping_tool/centromere_genotyping/samplesheet.csv --outdir outdir/ -resume

Draft 
Runs centromere analysis and outputs centromere architecture and mean aSat length per chromosome for short read data in different formats (raw reads, alignments or de-Brujin graphs). Works by counting the occurance of a predefined set of k-mers (using either jellyfish(LINK) for read level data or a simple grep based workflow for DBG data) and predicts centromeres by assessing presence and absence of k-mers in the predefined set and predict aSat length using k-mer dosage. For the length prediciton step the k-mer disage gets normalized using sets of k-mers unique (meaning n=1 per haplotype) across the human pangenome. Due to technical differences on how k-mers are found in raw short-read data vs. short-read assemblies (e.g. missing k-mers due to chance of truncated k-mers due to short-reads in short-read data; read correction might decrease number of k-mers while assembly should increase the number of k-mers due to contigs loonger than the inputted reads) there are differences depedning on the input data type. The prediciton works best with DBG transformed short-read data. To facilitate DBG generation the used workflow for DBG generation is located in this repo [LINK] and will be integrated as a optional step for raw reads in the nextflow pipeline soon. 

INPUT FORMAT
Output Format? 

ToDo
Add Input Columns for filtering DBG NORM?
Then perform filtering and add filtering information to output, dummy for short read
Column for DBG conversion of short read data
