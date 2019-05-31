import os
from pathlib import Path
THREADS = config["threads"]
TRIMMER = config["TRIMMER"]
ALIGNER = config["ALIGNER"]
METHOD = config["METHOD"]

rule all:
    input:
        expand("intermediate/genome.{n}.bt2", n = ["1","2","3","4"]),
        expand("intermediate/genome.rev.{n}.bt2", n = ["1","2"]),
        expand("intermediate/{sra_id}.{trimmer}.{aligner}.sam", sra_id=config["sample_ids"], trimmer=TRIMMER, aligner=ALIGNER),
        expand("intermediate/{sra_id}.{trimmer}.{aligner}.bam", sra_id=config["sample_ids"], trimmer=TRIMMER, aligner=ALIGNER),
        expand("intermediate/{sra_id}.{trimmer}.{aligner}.sorted.bam", sra_id=config["sample_ids"], trimmer=TRIMMER, aligner=ALIGNER)

rule bt2_index_genome:
    """
    Index a genome using Bowtie 2.
    """
    input:
        expand("genome/{genome_base}/{genome_base}_genomic.fna",
            genome_base = config["genomes"][config["genome_id"]]["url"].split("/")[-1])
    output:
        expand("intermediate/genome.{n}.bt2", n = ["1","2","3","4"]),
        expand("intermediate/genome.rev.{n}.bt2", n = ["1","2"])
    log:
        "log/index_genome.log"
    benchmark:
        "benchmarks/index_genome.bt2.index.benchmark.txt"
    shell:
        """
        bowtie2-build {input} intermediate/genome 2> {log}
        """

rule bt2_align:
    """
    Align a fastq file to a genome index using Bowtie 2.
    """
    input:
        fastq_1 = "transcriptome/reads/{sra_id}_1.{trimmer}.fastq",
        fastq_2 = "transcriptome/reads/{sra_id}_2.{trimmer}.fastq",
        index = expand("intermediate/genome.{n}.bt2", n = ["1","2","3","4"]),
        index_rev = expand("intermediate/genome.rev.{n}.bt2", n = ["1","2"])
    output:
        sam = temp("intermediate/{sra_id}.{trimmer}.bt2.sam"),
        bam = temp("intermediate/{sra_id}.{trimmer}.bt2.bam")
    threads: config["threads"]
    log:
        "log/{sra_id}.{trimmer}.bt2.align.log"
    benchmark:
        "benchmarks/{sra_id}.{trimmer}.bt2.align.benchmark.txt"
    run:
        # This gives the base name for the genome index, i.e. "intermediate/some_id"
        # rather than "intermediate/some_id.*.bt2"
        indexBase = "intermediate/genome"
        shell("bowtie2 -x " + indexBase + " --threads {threads} -1 {input.fastq_1} -2 {input.fastq_2} -S {output.sam} 2> {log}")
        shell("samtools view -bS {output.sam} > {output.bam}")

rule bbmap_align:
    params:
        genome_base = config["genomes"][config["genome_id"]]["url"].split("/")[-1]
    input:
        fastq_1 = "transcriptome/reads/{sra_id}_1.{trimmer}.fastq",
        fastq_2 = "transcriptome/reads/{sra_id}_2.{trimmer}.fastq",
        genome = expand("genome/{genome_base}/{genome_base}_genomic.fna", genome_base=config["genomes"][config["genome_id"]]["url"].split("/")[-1])
    output:
        sam = temp("intermediate/{sra_id}.{trimmer}.bbmap.sam"),
        bam = temp("intermediate/{sra_id}.{trimmer}.bbmap.bam")
    threads: THREADS
    log:
        "log/{sra_id}.{trimmer}.bbmap.align.log"
    benchmark:
        "benchmarks/{sra_id}.{trimmer}.bbmap.align.benchmark.txt"
    run:
        shell("bbmap.sh -Xmx20g trimreaddescriptions=t threads={threads} in1={input.fastq_1} in2={input.fastq_2} \
        out={output.sam} ref={input.genome}  path=genome/bbmap_index/ 2> {log}")
        shell("samtools view -bS {output.sam} > {output.bam}")


rule sort_bam:
    """
    Sort a bam file.
    """
    input:
        "intermediate/{sra_id}.{trimmer}.{aligner}.bam"
    output:
        "intermediate/{sra_id}.{trimmer}.{aligner}.sorted.bam"
    threads: THREADS
    shell:
        """
        samtools sort {input} > {output}
        """
