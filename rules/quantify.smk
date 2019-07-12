import os
from pathlib import Path
THREADS = config["THREADS"]
TRIMMER = config["TRIMMER"]
ALIGNER = config["ALIGNER"]
METHOD = config["METHOD"]

# Select genome/transcriptome ref_gff directories
if genome_url:
    ref_gff = expand("reference/genome/{genome_id}/{genome_id}_genomic.gff",
        genome_id = config["genomes"][config["genome_id"]]["url"].split("/")[-1])
    ref_fna = expand("reference/genome/{genome_id}/{genome_id}_cds_from_genomic_salmon.fna",
        genome_id = config["genomes"][config["genome_id"]]["url"].split("/")[-1])
else:
    if ASSEMBLER=="megahit":
        ref_gff = "reference/assembled/megahit_out/genes_annotated.gff"
        ref_fna = "reference/assembled/megahit_out/genes_annotated.fna"
    elif ASSEMBLER=="trinity":
        ref_gff = "reference/assembled/trinity_out/genes_annotated.gff"
        ref_fna = "reference/assembled/trinity_out/genes_annotated.fna"

rule htseq_count_table:
    """
    Generate a count table using htseq-count.
    """
    input:
        bams = expand("intermediate/{sra_id}.{trimmer}.{aligner}.sorted.bam",
        sra_id = config["sample_ids"], trimmer = config["TRIMMER"], aligner = config["ALIGNER"]),
        gff = ref_gff
    output:
        "results/tables/htseq.{aligner}.{trimmer}.counts.tsv"
    shadow: "minimal"
    threads: THREADS
    benchmark:
        "benchmarks/htseq.{aligner}.{trimmer}.benchmark.txt"
    shell:
        """
        # Save the count table as a temporary file and then prepend a header line
        # with the sample names
        htseq-count --format bam --type gene --idattr locus_tag -r pos {input.bams} {input.gff} > {output}
        #echo '{input.bams}' | tr ' ' '\t' | cat - tempfile2 > {output}
        """

rule feature_counts_table:
    """
    Generate a count table with featureCounts
    """
    input:
        bams=expand("intermediate/{sra_id}.{trimmer}.{aligner}.sorted.bam",
        sra_id = config["sample_ids"], trimmer = config["TRIMMER"], aligner = config["ALIGNER"]),
        gff = ref_gff
    output:
        "results/tables/featureCounts.{aligner}.{trimmer}.counts.tsv"
    threads: THREADS
    log:
        "log/featureCounts.{aligner}.{trimmer}.log"
    benchmark:
        "benchmarks/featureCounts.{aligner}.{trimmer}.benchmark.txt"
    shell:
        """
        featureCounts -O -p -T {threads} -t gene -F GTF -a {input.gff} -o {output} {input.bams}
        """

rule salmon_index:
    input:
        fasta = ref_fna
    output:
        directory("reference/salmon_quasi")
    threads: THREADS
    run:
        shell("salmon index -t {input.fasta} -i {output} --type quasi -k 31")


rule salmon_quant:
    """
    Generate directories containing count files with salmon (quasi mode)
    """
    input:
        fastq_1="transcriptome/reads/{sra_id}_1.{trimmer}.fastq",
        fastq_2="transcriptome/reads/{sra_id}_2.{trimmer}.fastq",
        salmon_dir = directory("reference/salmon_quasi")
    output:
        directory("transcriptome/salmon/{sra_id}_{trimmer}")
    threads: THREADS
    log:
        "log/salmon.{sra_id}.{trimmer}.log"
    benchmark:
        "benchmarks/salmon.quant.{sra_id}.{trimmer}.benchmark.txt"
    run:
        shell("salmon quant -i {input.salmon_dir} -l A -p {threads} --validateMappings \
        -1 {input.fastq_1} -2 {input.fastq_2} -o {output} 2> {log}")


rule salmon_quant_table:
    input:
        expand(directory("transcriptome/salmon/{sra_id}_{trimmer}"), sra_id = config["sample_ids"], trimmer = TRIMMER)
    output:
        "results/tables/salmon.{trimmer}.counts.tsv"
    run:
        shell("salmon quantmerge --quants {input} --column numreads -o {output}")
