# PRAXIS: PRokaryotic Activity and eXpression Informatics System
End-to-end de novo transcriptomic pipeline
v0.2

Authors: David Levy-Booth and Parker Lloyd


This snakemake pipeline aids in downloading and processing of genomic and transcriptomic sequences, transcriptomic assembly and quantification using modularized, benchmarked tools. 

Steps to Initialize:

Download or clone repository:

    git clone https://github.com/davidlevybooth/Praxis.git

Change to the Praxis directory:

    cd Praxis

Create the conda environment:

    conda env create --name Praxis --file envs/environment.yaml
    
Activate conda environment:

    conda activate Praxis

Install Praxis:

    python setup.py install
    
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

    -j    --jobs        The number of jobs. Analogous to snakemake -j#
 
    -r     --read_dir        The directory in which fastq reads are stored
    
    -o    --outdir        The directory for the final outputs
    
    -c    --configfile         The yaml file containing all of the above configuration details
    

Declaring a reference:

There are 3 ways to provide a reference to the pipeline. The path to a genome file can be provided using the corresponding flag (-f), the NCBI URL can be provided for the pipeline to download the genome (-u), or the pipeline can use the provided RNASeq reads located in the directory transcriptome/reads/untrimmed to assemble a reference transcriptome.

If the path to a user provided genome file is specified, it will always be used as the reference, regardless of any other options available. To download a reference from NCBI, remove the path to the genome file and specify the NCBI URL (Praxis -f -u [ncbi_url]). To assemble a reference transcriptome, make sure both the genome file and NCBI URL are not specified, make sure there are RNASeq reads located in the directory specified in transcriptome/reads/untrimmed, and make sure the desired assembler is chosen.
 
User Provided RNASeq Reads:

If you are providing your own RNASeq fastq files, the filenames will need a certain format. The pipeline expects a particular naming convention for forward and reverse RNASeq reads. Forward reads must have the form {xxx}_1.fastq and reverse reads must have the form {xxx}_2.fastq. 

The columns in the samples.tsv file must be updated to include the file names for each of the forward and reverse reads. 
 
 
