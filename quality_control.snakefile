import os
from pathlib import Path
THREADS = config["threads"]
TRIMMER = config["TRIMMER"]
ALIGNER = config["ALIGNER"]
METHOD = config["METHOD"]
if "salmon" in METHOD and len(METHOD) == 1:
    count_out = "results/tables/salmon.{{trimmer}}.counts.tsv"
if "salmon" in METHOD and len(METHOD) > 1:
    METHOD.remove("salmon")
    count_out = ["results/tables/salmon.{trimmer}.counts.tsv",
    "results/tables/{method}.{aligner}.{trimmer}.counts.tsv"]

rule all:
    input:
       expand("transcriptome/reads/{sra_id}_1.trimmomatic.fastq", sra_id = config["sample_ids"]),
       expand("transcriptome/reads/{sra_id}_2.trimmomatic.fastq", sra_id = config["sample_ids"]),
       expand("transcriptome/reads/{sra_id}_1.unpaired.fastq", sra_id = config["sample_ids"]),
       expand("transcriptome/reads/{sra_id}_2.unpaired.fastq", sra_id = config["sample_ids"]),
       expand("transcriptome/reads/{sra_id}_1.bbduk.fastq", sra_id = config["sample_ids"]),
       expand("transcriptome/reads/{sra_id}_2.bbduk.fastq", sra_id = config["sample_ids"]),
       expand("transcriptome/qc/fastqc/{trimmer}/{sra_id}_1.html", sra_id = config["sample_ids"], trimmer = TRIMMER),
       expand("transcriptome/qc/fastqc/{trimmer}/{sra_id}_2.html", sra_id = config["sample_ids"], trimmer = TRIMMER),
       expand(directory("transcriptome/qc/fastqc/{trimmer}/{sra_id}_1"), sra_id = config["sample_ids"], trimmer = TRIMMER),
       expand(directory("transcriptome/qc/fastqc/{trimmer}/{sra_id}_2"), sra_id = config["sample_ids"], trimmer = TRIMMER)


rule trimmomatic:
    input:
        read1 = "transcriptome/reads/{sra_id}_1.untrimmed.fastq",
        read2 = "transcriptome/reads/{sra_id}_2.untrimmed.fastq"
    output:
       paired1 = "transcriptome/reads/{sra_id}_1.trimmomatic.fastq",
       paired2 = "transcriptome/reads/{sra_id}_2.trimmomatic.fastq",
       unpaired1 = "transcriptome/reads/{sra_id}_1.unpaired.fastq",
       unpaired2 = "transcriptome/reads/{sra_id}_2.unpaired.fastq"
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
        read1 = "transcriptome/reads/{sra_id}_1.untrimmed.fastq",
        read2 = "transcriptome/reads/{sra_id}_2.untrimmed.fastq"
    output:
       paired1 = "transcriptome/reads/{sra_id}_1.bbduk.fastq",
       paired2 = "transcriptome/reads/{sra_id}_2.bbduk.fastq"
    threads: THREADS
    log:
        "log/{sra_id}.bbduk.log"
    benchmark:
        "benchmarks/{sra_id}.bbduk.benchmark.txt"
    shell:
        """
        bbduk.sh -Xmx20g t={threads} in1={input.read1} in2={input.read2} out1={output.paired1} out2={output.paired2} \
        ref=utils/adapters/TruSeq3-PE-2.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=15 2> {log}
        """

rule fastqc:
    input:
        "transcriptome/reads/{sra_id}_{num}.{trimmer}.fastq"
    output:
        html="transcriptome/qc/fastqc/{trimmer}/{sra_id}_{num}.html",
        unzipped=directory("transcriptome/qc/fastqc/{trimmer}/{sra_id}_{num}")
    params: ""
    wildcard_constraints:
        num="[1,2]"
    log:
        "log/fastqc/{trimmer}/{sra_id}_{num}.log"
    script:
        "fastqc.py"
