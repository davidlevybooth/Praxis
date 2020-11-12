#!/usr/bin/env python

__authors__ = "David Levy-Booth, Parker Lloyd"
__copyright__ = "Copyright 2019, David Levy-Booth"
__email__ = "dlevyboo@mail.ubc.ca"
__license__ = "GPL3"
#!/usr/bin/python3

import sys
import os
from os import listdir
from os.path import isfile, join
import itertools
import pandas as pd
from pathlib import Path
import subprocess
import yaml
import logging
import argparse
import psutil
import copy
import time
from datetime import datetime
now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")

print("\033[94m\n  ╔═══════════════════════════════════════════════════════════╗\033[0m")
print("\033[94m  ║      _____    ______  _______  _     _  _____  _______    ║\033[0m")
print("\033[94m  ║     |_____]  |_____/  |_____|   \\___/     |    |______    ║\033[0m")
print("\033[94m  ║     |        |    \\_  |     |  _/   \\_  __|__  ______|    ║\033[0m")
print("\033[94m  ║                                                           ║\033[0m")
print("\033[94m  ║ Prokaryotic Activity and Expression Informatics System    ║\033[0m")
print("\033[94m  ║ v0.2 | D. Levy-Booth | UBC Sequencing and Bioinformatics  ║\033[0m")
print("\033[94m  ║                                                           ║\033[0m")
print("\033[94m  ║ GitHub:  https://github.com/davidlevybooth/Praxis         ║\033[0m")
print("\033[94m  ║ Twitter: https://twitter.com/davidlevybooth               ║\033[0m")
print("\033[94m  ║                                                           ║\033[0m")
print("\033[94m  ║ `Praxis -h` for help                                      ║\033[0m")
print("\033[94m  ║                                                           ║\033[0m")
print("\033[94m  ╚═══════════════════════════════════════════════════════════╝\n\033[0m")
time.sleep(1)

def get_snakefile(file="Snakefile"):
    sf = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), file)
    if not os.path.exists(sf):
        sys.exit("Unable to locate the Snakemake workflow file; tried %s" % sf)
    return sf


def main():
    # Define flags
    parser = argparse.ArgumentParser(description='pipeline configuration')
    parser.add_argument("-d", "--method", default = "salmon", type = str, help = "Either featureCounts, htseq, or salmon.")
    parser.add_argument("-a", "--aligner", default = "none", type = str, help = "Either bt2, bbmap, or none.")
    parser.add_argument("-q", "--trimmer", default = "trimmomatic", type = str, help = "Either trimmomatic, bbduk, or untrimmed.")
    parser.add_argument("-s", "--assembler", default = "megahit", type = str, help = "Either trinity or megahit.")
    parser.add_argument("-t", "--threads", default = 2, type = int, help = "The number of threads the pipeline is allowed to use.")
    parser.add_argument("-m", "--max_memory", type = int, help = "The amount of memory provided to the assembler. Enter in byte format.")
    parser.add_argument("-u", "--genome_url", type = str, help = "The NCBI url from which the reference genome can be downloaded.")
    parser.add_argument("-f", "--genome_file", type = str, help = "The file containing the reference genome.")
    parser.add_argument("-r", "--read_dir", type = str, help = "The directory in which fastq reads are stored")
    parser.add_argument("-o", "--outdir", type = str, help = "The directory for the final outputs.")
    parser.add_argument("-c", "--configfile", default = "config.yaml", type = str, help = "yaml file containing all of the above configuration details.")
    parser.add_argument("-j", "--jobs", default = 2, type = str, help = "number of jobs")
    args = parser.parse_args()


    #Assign arguments 
    # Iterate through all possible flags and update corresponding field in config image if flag is present and value is accepted
    if args.configfile:
        configfile = args.configfile
    if args.jobs:
        jobs = args.jobs

    # Open and create an image of the Snakemake config
    stream = open(configfile, "r")
    original_data = yaml.load(stream, yaml.FullLoader)
    new_data = copy.deepcopy(original_data)

    # else:
    for arg in vars(args):
        if getattr(args, arg) is not None:
            if arg == "method":
                val = getattr(args, arg)
                if val in ["featureCounts", "htseq", "salmon"]:
                    METHOD = val
                else:
                    raise Exception("Requested method is not an available option.")
            elif arg == "aligner":
                val = getattr(args, arg)
                if val in ["bt2", "bbmap", "none"]:
                    ALIGNER = val
                else:
                    raise Exception("Requested aligner is not an available option.")
            elif arg == "trimmer":
                val = getattr(args, arg)
                if val in ["trimmomatic", "bbduk", "untrimmed"]:
                    TRIMMER = val
                else:
                    raise Exception("Requested trimmer is not an available option.")
            elif arg == "assembler":
                val = getattr(args, arg)
                if val in ["trinity", "megahit", "none"]:
                    ASSEMBLER = val
                else:
                    raise Exception("Requested assembler is not an available option.")
            elif arg == "threads":
                val = getattr(args, arg)
                if val <= psutil.cpu_count():
                    THREADS = val
                else:
                    raise Exception("Number of threads requested exceeds number of available cores.")
            elif arg == "max_memory":
                val = getattr(args, arg)
                if val <= psutil.virtual_memory()[1]:
                    MAX_MEMORY = val
                else:
                    raise Exception("Requested memory exceeds available memory")
           
            elif arg == "genome_url":
                genome_url = getattr(args, arg)
            elif arg == "genome_file":
                if getattr(args, arg).split(".")[-1] in ["fasta", "fa", "fna"]:
                    genome_file = getattr(args, arg)
                else:
                    raise Exception("Provided reference file is not in an accepted format (.fasta, .fa, .fna).")

    # Rewrite config file if changes were made to image
    if original_data != new_data:
        with open(configfile, 'w') as yaml_file:
            yaml_file.write(yaml.dump(new_data, default_flow_style=False))

    # Close file stream
    stream.close()

    #Actually running snakemake
    cmd = (
        "snakemake -s {snakefile} "
        "-j{jobs} --rerun-incomplete "
        "--configfile {configfile}"
    ).format(
        snakefile=get_snakefile(),
        jobs=jobs,
        configfile=configfile),
    logging.info("Executing: %s" % cmd)

    print(cmd)

    try:    
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        # removes the traceback
        logging.critical(e)
        exit(1)

if __name__ == ' __main__':
    main()