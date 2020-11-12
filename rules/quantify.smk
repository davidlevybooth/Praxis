import os
from pathlib import Path
THREADS = config["THREADS"]
TRIMMER = config["TRIMMER"]
ALIGNER = config["ALIGNER"]
METHOD = config["METHOD"]

genome_url = config["genome"]["ncbi_url"]
genome_file = config["genome"]["ref_file"]
if genome_file:
    ref_dir = "/".join(genome_file.split("/")[:-2])

# Select genome/transcriptome ref_gff directories

#Determine if annotation required
if genome_file:
    feature_type = "gene"
    if config["genome"]["annotate"]:
        ref_gff = ref_dir + "/genes_annotated.gff"
        ref_fna = ref_dir + "/genes_annotated.fna"
    else: 
        ref_fna = genome_file
        ref_gff = genome_file.replace(".fna", ".gff")

#If annotation, set file paths
else:
    if genome_url:
        feature_type = "gene"
        ref_gff = expand("reference/genome/{genome_id}/{genome_id}_genomic.gff",
            genome_id = config["genome"]["ncbi_url"].split("/")[-1])
        ref_fna = expand("reference/genome/{genome_id}/{genome_id}_cds_from_genomic_salmon.fna",
            genome_id = config["genome"]["ncbi_url"].split("/")[-1])
    else:
        feature_type = "CDS"
        if ASSEMBLER=="megahit":
            ref_gff = "reference/assembled/megahit_genes/genes_annotated.gff"
            ref_fna = "reference/assembled/megahit_genes/genes_annotated.fna"
        elif ASSEMBLER=="trinity":
            ref_gff = "reference/assembled/trinity_genes/genes_annotated.gff"
            ref_fna = "reference/assembled/trinity_genes/genes_annotated.fna"

rule htseq_count_table:
    """
    Generate a count table using htseq-count.
    """
    input:
        bams = expand(output_directory + "intermediate/{trimmer}/{aligner}/{sra_id}.sorted.bam", trimmer = config["TRIMMER"], aligner = config["ALIGNER"], sra_id = sampleID_list),
        gff = ref_gff
    output:
        output_directory + "results/tables/htseq/{trimmer}/{aligner}/counts.tsv"
    shadow: "minimal"
    threads: THREADS
    benchmark:
        output_directory + "benchmarks/htseq/{trimmer}/{aligner}/benchmark.txt"
    shell:
        """
        # Save the count table as a temporary file and then prepend a header line
        # with the sample names
        htseq-count --format bam --type gene --idattr gene_id -r pos {input.bams} {input.gff} > {output}
        #echo '{input.bams}' | tr ' ' '\t' | cat - tempfile2 > {output}
        """

rule feature_counts_table:
    """
    Generate a count table with featureCounts
    """
    input:
        bams = expand(output_directory + "intermediate/{trimmer}/{aligner}/{sra_id}.sorted.bam", trimmer=TRIMMER, aligner=ALIGNER, sra_id = sampleID_list),
        gff = ref_gff
    output:
        output_directory + "results/tables/featureCounts/{trimmer}/{aligner}/counts.tsv"
    params:
        type = feature_type
    threads: THREADS
    log:
       output_directory + "log/featureCounts/{trimmer}/{aligner}/log"
    benchmark:
        output_directory + "benchmarks/featureCounts/{trimmer}/{aligner}/benchmark.txt"
    shell:
        """
        featureCounts -O -p -T {threads} -t {params.type} -F GTF -a {input.gff} -o {output} {input.bams}
        """

rule salmon_index:
    """
    Index a reference with salmon.
    """
    input:
        fasta = ref_fna
    output:
        directory(output_directory + "reference/salmon_quasi"),
    threads: THREADS
    run:
        shell("salmon index -t {input.fasta} -i {output} --type quasi -k 31")


rule salmon_quant:
    """
    Generate directories containing count files with salmon (quasi mode).
    """
    input:
        fastq_1= output_directory + "transcriptome/reads/{trimmer}/{sra_id}_1.fastq",
        fastq_2=output_directory + "transcriptome/reads/{trimmer}/{sra_id}_2.fastq",
        salmon_dir = output_directory + "reference/salmon_quasi"
    output:
        directory(output_directory + "transcriptome/salmon/{sra_id}_{trimmer}")
    threads: THREADS
    log:
        "log/salmon/{trimmer}/{sra_id}.log"
    benchmark:
        "benchmarks/salmon.quant/{trimmer}/{sra_id}.benchmark.txt"
    run:
        shell("salmon quant -i {input.salmon_dir} -l A -p {threads} --validateMappings \
        -1 {input.fastq_1} -2 {input.fastq_2} -o {output} 2> {log}")


rule salmon_quant_table:
    """
    Generate a count table with salmon.
    """
    input:
        expand(output_directory + "transcriptome/salmon/{sra_id}_{trimmer}", trimmer = TRIMMER, sra_id = sampleID_list)
    output:
        output_directory + "results/tables/salmon/{trimmer}/counts.tsv"
    run:
        shell("salmon quantmerge --quants {input} --column numreads -o {output}")
