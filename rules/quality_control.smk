import os
import pandas as pd
from pathlib import Path

# Verify that the RNASeq/SRA files exist

if ("sample_ids" in config):
    sample_ids = config["sample_ids"]
else: 
    reads = pd.read_csv('samples.tsv', delimiter = '\t')['Forward_Reads'].tolist() + pd.read_csv('samples.tsv', delimiter = '\t')['Reverse_Reads'].tolist()
    #path = os.getcwd() + sample_directory
    path = sample_directory
    for file in reads:
        if not os.path.isfile(path + file):
            raise Exception('\'' + file + '\' is not present in \'' + path + '\'.')

rule trimmomatic:
    """
    Trim RNASeq fasta file using trimmomatic.
    """
    input:
        # read1 = sample_directory + "{sra_id}_1.fastq",
        # read2 = sample_directory + "{sra_id}_2.fastq"
        read1 = lambda wildcards: sample_directory + samples_table.Forward_Reads[samples_table.SampleID == wildcards.sampleID],
        read2 = lambda wildcards: sample_directory + samples_table.Reverse_Reads[samples_table.SampleID == wildcards.sampleID],
    output:
       paired1 = output_directory + "transcriptome/reads/trimmomatic/{sampleID}_1.fastq",
       paired2 = output_directory + "transcriptome/reads/trimmomatic/{sampleID}_2.fastq",
       unpaired1 = output_directory + "transcriptome/reads/unpaired/{sampleID}_1.fastq",
       unpaired2 = output_directory + "transcriptome/reads/unpaired/{sampleID}_2.fastq"
    params:
        adapter = config["adapter_file"]
    threads: THREADS
    log:
        output_directory + "log/{sampleID}.trimmomatic.log"
    benchmark:
        output_directory + "benchmarks/{sampleID}.trimmomatic.benchmark.txt"
    shell:
      "trimmomatic PE -threads {threads} -phred33 "
      "{input.read1} {input.read2} {output.paired1} {output.unpaired1} {output.paired2} {output.unpaired2} "
      "ILLUMINACLIP:{params.adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2> {log}"

rule bbduk:
    """
    Trim RNASeq fasta file using bbduk.
    """
    input:
        read1 = sample_directory + "{sra_id}_1.fastq",
        read2 = sample_directory + "{sra_id}_2.fastq"
    output:
       paired1 = output_directory + "transcriptome/reads/bbduk/{sra_id}_1.fastq",
       paired2 = output_directory + "transcriptome/reads/bbduk/{sra_id}_2.fastq"
    params:
        adapter = config["adapter_file"]
    threads: THREADS
    log:
        output_directory + "log/bbduk/{sra_id}.log"
    benchmark:
        output_directory + "benchmarks/bbduk/{sra_id}.benchmark.txt"
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
        output_directory + "transcriptome/reads/trimmed/" + "{sra_id}_{num}.fastq"
    output:
        html= output_directory + "transcriptome/qc/fastqc/{trimmer}/{sra_id}_{num}.html",
        unzipped=directory(output_directory + "transcriptome/qc/fastqc/{trimmer}/{sra_id}_{num}")
    wildcard_constraints:
        num="[1,2]"
    log:
        output_directory + "log/fastqc/{trimmer}/{sra_id}_{num}.log"
    script:
        "../scripts/fastqc.py"
