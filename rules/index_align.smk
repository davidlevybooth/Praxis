import os
from pathlib import Path
THREADS = config["THREADS"]
TRIMMER = config["TRIMMER"]
ALIGNER = config["ALIGNER"]
METHOD = config["METHOD"]
ASSEMBLER = config["ASSEMBLER"]

genome_url = config["genomes"][config["genome_id"]]["url"]

# Select genome/transcriptome reference directories
if genome_url:
    reference = expand("reference/genome/{genome_id}/{genome_id}_genomic.fna",
        genome_id = genome_url.split("/")[-1])
    indexBase = "intermediate/genome"
else:
    if ASSEMBLER=="megahit":
        reference = "reference/assembled/megahit_out/final.contigs.fa"
        indexBase = "intermediate/megahit/transcriptome"
    elif ASSEMBLER=="trinity":
        reference = "reference/assembled/trinity_out/Trinity.fasta"
        indexBase = "intermediate/trinity/transcriptome"

if "bt2" in ALIGNER:
    rule bt2_index:
        """
        Index a genome using Bowtie 2.
        """
        input:
            ref=reference
        output:
            expand("{indexBase}.{n}.bt2", indexBase = indexBase, n = ["1","2","3","4"]),
            expand("{indexBase}.rev.{n}.bt2", indexBase = indexBase, n = ["1","2"])
        params:
            base = indexBase
        log:
            "log/index.log"
        benchmark:
            "benchmarks/index.bt2.index.benchmark.txt"
        shell:
            """
            bowtie2-build {input} {params.base} 2> {log}
            """

    rule bt2_align:
        """
        Align a fastq file to a genome index using Bowtie 2.
        """
        input:
            fastq_1 = "transcriptome/reads/{sra_id}_1.{trimmer}.fastq",
            fastq_2 = "transcriptome/reads/{sra_id}_2.{trimmer}.fastq",
            index = expand("{indexBase}.{n}.bt2", indexBase = indexBase, n = ["1","2","3","4"]),
            index_rev = expand("{indexBase}.rev.{n}.bt2", indexBase = indexBase, n = ["1","2"])
        output:
            sam = temp("intermediate/{sra_id}.{trimmer}.bt2.sam"),
            bam = temp("intermediate/{sra_id}.{trimmer}.bt2.bam")
        params:
            base = indexBase
        threads: THREADS
        log:
            "log/{sra_id}.{trimmer}.bt2.align.log"
        benchmark:
            "benchmarks/{sra_id}.{trimmer}.bt2.align.benchmark.txt"
        run:
            shell("bowtie2 -x " + params.base + " --threads {threads} -1 {input.fastq_1} -2 {input.fastq_2} -S {output.sam} 2> {log}")
            shell("samtools view -bS {output.sam} > {output.bam}")

if "bbmap" in ALIGNER:
    rule bbmap_align:
        input:
            fastq_1 = "transcriptome/reads/{sra_id}_1.{trimmer}.fastq",
            fastq_2 = "transcriptome/reads/{sra_id}_2.{trimmer}.fastq",
            ref = reference
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
            out={output.sam} ref={input.ref}  path=reference/bbmap_index/ 2> {log}")
            shell("samtools view -bS {output.sam} > {output.bam}")

rule sort_bam:
    """
    Sort a bam file.
    """
    input:
        "intermediate/{sra_id}.{trimmer}.{aligner}.bam"
    output:
        temp("intermediate/{sra_id}.{trimmer}.{aligner}.sorted.bam")
    threads: THREADS
    shell:
        """
        samtools sort {input} > {output}
        """
