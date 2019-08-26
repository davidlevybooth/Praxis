import os
import pandas as pd
from pathlib import Path
THREADS = config["THREADS"]

# Verify that the RNASeq/SRA files exist
if not config["sample_ids"]:
    reads = pd.read_csv('samples.tsv', delimiter = '\t')['Forward_Reads'].tolist() + pd.read_csv('samples.tsv', delimiter = '\t')['Reverse_Reads'].tolist()
    path = os.getcwd() + '/transcriptome/reads/untrimmed/'
    for file in reads:
        if not os.path.isfile(path + file):
            raise Exception('\'' + file + '\' is not present in \'' + path + '\'.')

rule trimmomatic:
    """
    Trim RNASeq fasta file using trimmomatic.
    """
    input:
        read1 = "transcriptome/reads/untrimmed/{sra_id}_1.fastq",
        read2 = "transcriptome/reads/untrimmed/{sra_id}_2.fastq"
    output:
       paired1 = "transcriptome/reads/trimmomatic/{sra_id}_1.fastq",
       paired2 = "transcriptome/reads/trimmomatic/{sra_id}_2.fastq",
       unpaired1 = "transcriptome/reads/unpaired/{sra_id}_1.fastq",
       unpaired2 = "transcriptome/reads/unpaired/{sra_id}_2.fastq"
    params:
        adapter = config["adapter_file"]
    threads: THREADS
    log:
        "log/{sra_id}.trimmomatic.log"
    benchmark:
        "benchmarks/{sra_id}.trimmomatic.benchmark.txt"
    shell:
      "trimmomatic PE -threads {threads} -phred33 "
      "{input.read1} {input.read2} {output.paired1} {output.unpaired1} {output.paired2} {output.unpaired2} "
      "ILLUMINACLIP:{params.adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2> {log}"

rule bbduk:
    """
    Trim RNASeq fasta file using bbduk.
    """
    input:
        read1 = "transcriptome/reads/untrimmed/{sra_id}_1.fastq",
        read2 = "transcriptome/reads/untrimmed/{sra_id}_2.fastq"
    output:
       paired1 = "transcriptome/reads/bbduk/{sra_id}_1.fastq",
       paired2 = "transcriptome/reads/bbduk/{sra_id}_2.fastq"
    params:
        adapter = config["adapter_file"]
    threads: THREADS
    log:
        "log/bbduk/{sra_id}.log"
    benchmark:
        "benchmarks/bbduk/{sra_id}.benchmark.txt"
    shell:
        """
        bbduk.sh -Xmx20g t={threads} in1={input.read1} in2={input.read2} out1={output.paired1} out2={output.paired2} \
        ref={params.adapter} ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=15 2> {log}
        """

rule fastqc:
    """
    Store quality scores for sequenced reads.
    """
    input:
        "transcriptome/reads/{trimmer}/{sra_id}_{num}.fastq"
    output:
        html="transcriptome/qc/fastqc/{trimmer}/{sra_id}_{num}.html",
        unzipped=directory("transcriptome/qc/fastqc/{trimmer}/{sra_id}_{num}")
    wildcard_constraints:
        num="[1,2]"
    log:
        "log/fastqc/{trimmer}/{sra_id}_{num}.log"
    script:
        "scripts/fastqc.py"
