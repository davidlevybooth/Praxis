# PRAXIS: PRokaryotic Activity and eXpression Informatics System
End-to-end de novo transcriptomic pipeline

Authors: David Levy-Booth and Parker Lloyd


This snakemake pipeline aids in downloading and processing of genomic and transcriptomic sequences, transcriptomic assembly and quantification using modularized, benchmarked tools. 

Steps to Initialize:

Download or clone repository:

git clone https://github.com/davidlevybooth/Praxis.git

Run init.sh to create Praxis function and conda environment:

    ./init.sh
    
Activate conda environment:

    conda activate <env_name>
    
Run pipeline with desired options:

    Praxis [ -d featureCounts -a bbmap -q bbduk -s trinity … ]
    
 
Pipeline Options:

    -d    --method    The tool that will be used to generate the count table of expressed genes. The available tools are featureCounts, htseq, and salmon
    
    -a    --aligner    The tool that will be used to align the RNASeq reads to the reference genome/transcriptome
    
    -q    --trimmer    The tool that will trim contaminant sequences and adapter sequences off the RNASeq reads
    
    -s    --assembler    The tool that will assemble a reference transcriptome if there is no genome URL provided
    
    -t    --threads    The number of threads that the pipeline is allowed to use.
    
    -m    --max_memory    The amount of memory provided to the assembler. Enter in byte format
    
    -u    --genome_url        The NCBI url from which the reference genome can be downloaded
    
    -f    --genome_file        The file containing the reference genome
 
Declaring a reference

There are 3 ways to provide a reference to the pipeline. The path to a genome file can be provided using the corresponding flag (-f), the NCBI URL can be provided for the pipeline to download the genome (-u), or the pipeline can use the provided RNASeq reads located in the directory transcriptome/reads/untrimmed to assemble a reference transcriptome.

If the path to a user provided genome file is specified, it will always be used as the reference, regardless of any other options available. To download a reference from NCBI, remove the path to the genome file and specify the NCBI URL (Praxis -f -u [ncbi_url]). To assemble a reference transcriptome, make sure both the genome file and NCBI URL are not specified, make sure there are RNASeq reads located in the directory specified in transcriptome/reads/untrimmed, and make sure the desired assembler is chosen.
 
User Provided RNASeq Reads

If you are providing your own RNASeq fastq files, the filenames will need a certain format. The pipeline expects a particular naming convention for forward and reverse RNASeq reads. Forward reads must have the form {xxx}_1.fastq and reverse reads must have the form {xxx}_2.fastq. 

Each of the different strings {xxx} must be listed in the config parameter sample_ids. 

These RNASeq fastq files must be placed in a directory called transcriptome/reads/untrimmed. Create this directory if it does not already exist.

The columns in the samples.tsv file must be updated to include the file names for each of the forward and reverse reads.
 
For future reference: Look into creating user defined output folders to organize multiple datasets
 
 
 
