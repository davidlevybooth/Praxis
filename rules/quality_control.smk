import os
from pathlib import Path
THREADS = config["THREADS"]
TRIMMER = config["TRIMMER"]
ALIGNER = config["ALIGNER"]
METHOD = config["METHOD"]

rule trimmomatic:
    input:
        read1 = "transcriptome/reads/untrimmed/{sra_id}_1.fastq",
        read2 = "transcriptome/reads/untrimmed/{sra_id}_2.fastq"
    output:
       paired1 = "transcriptome/reads/trimmomatic/{sra_id}_1.fastq",
       paired2 = "transcriptome/reads/trimmomatic/{sra_id}_2.fastq",
       unpaired1 = "transcriptome/reads/unpaired/{sra_id}_1.fastq",
       unpaired2 = "transcriptome/reads/unpaired/{sra_id}_2.fastq"
    threads: THREADS
    log:
        "log/{sra_id}.trimmomatic.log"
    benchmark:
        "benchmarks/{sra_id}.trimmomatic.benchmark.txt"
    shell:
      "trimmomatic PE -threads {threads} -phred33 "
      "{input.read1} {input.read2} {output.paired1} {output.unpaired1} {output.paired2} {output.unpaired2} "
      "ILLUMINACLIP:utils/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2> {log}"

rule bbduk:
    input:
        read1 = "transcriptome/reads/untrimmed/{sra_id}_1.fastq",
        read2 = "transcriptome/reads/untrimmed/{sra_id}_2.fastq"
    output:
       paired1 = "transcriptome/reads/bbduk/{sra_id}_1.fastq",
       paired2 = "transcriptome/reads/bbduk/{sra_id}_2.fastq"
    threads: THREADS
    log:
        "log/bbduk/{sra_id}.log"
    benchmark:
        "benchmarks/bbduk/{sra_id}.benchmark.txt"
    shell:
        """
        bbduk.sh -Xmx20g t={threads} in1={input.read1} in2={input.read2} out1={output.paired1} out2={output.paired2} \
        ref=utils/adapters/TruSeq3-PE-2.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=15 2> {log}
        """

rule fastqc:
    input:
        "transcriptome/reads/{trimmer}/{sra_id}_{num}.fastq"
    output:
        html="transcriptome/qc/fastqc/{trimmer}/{sra_id}_{num}.html",
        unzipped=directory("transcriptome/qc/fastqc/{trimmer}/{sra_id}_{num}")
    params: ""
    wildcard_constraints:
        num="[1,2]"
    log:
        "log/fastqc/{trimmer}/{sra_id}_{num}.log"
    script:
        "scripts/fastqc.py"
