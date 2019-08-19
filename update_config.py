import os
import yaml
import argparse
import psutil
from os import listdir
from os.path import isfile, join

parser = argparse.ArgumentParser(description='pipeline configuration')
parser.add_argument("-d", "--method", type = str, help = "Either featureCounts, htseq, or salmon.")
parser.add_argument("-a", "--aligner", type = str, help = "Either bt2, bbmap, or none.")
parser.add_argument("-q", "--trimmer", type = str, help = "Either trimmomatic, bbduk, or untrimmed.")
parser.add_argument("-s", "--assembler", type = str, help = "Either trinity or megahit.")
parser.add_argument("-t", "--threads", type = int, help = "The number of threads the pipeline is allowed to use.")
parser.add_argument("-m", "--max_memory", type = int, help = "The amount of memory provided to the assembler. Enter in byte format.")
parser.add_argument("-u", "--genome_url", type = str, help = "The NCBI url from which the reference genome can be downloaded.")
parser.add_argument("-f", "--genome_file", type = str, help = "The file containing the reference genome.")

args = parser.parse_args()

for arg in vars(args):
    # print(arg, getattr(args, arg))
    if getattr(args, arg) is not None:

        # open and read the Snakemake config
        file = "config.yaml"
        stream = open(file, "r")
        data = yaml.load(stream, yaml.FullLoader)

        if arg is "method":
            val = getattr(args, arg)
            if val in ["featureCounts", "htseq", "salmon"]:
                data[arg.upper()] = val
            else:
                raise Exception("Requested method is not an available option.")
        elif arg is "aligner":
            val = getattr(args, arg)
            if val in ["bt2", "bbmap", "none"]:
                data[arg.upper()] = val
            else:
                raise Exception("Requested aligner is not an available option.")
        elif arg is "trimmer":
            val = getattr(args, arg)
            if val in ["trimmomatic", "bbduk", "untrimmed"]:
                data[arg.upper()] = val
            else:
                raise Exception("Requested trimmer is not an available option.")
        elif arg is "assembler":
            val = getattr(args, arg)
            if val in ["trinity", "megahit"]:
                data[arg.upper()] = val
            else:
                raise Exception("Requested assembler is not an available option.")
        elif arg is "threads":
            val = getattr(args, arg)
            if val <= psutil.cpu_count():
                data[arg.upper()] = val
            else:
                raise Exception("Number of threads requested exceeds number of available cores.")
        elif arg is "max_memory":
            val = getattr(args, arg)
            if val <= psutil.virtual_memory()[1]:
                data[arg.upper()] = val
            else:
                raise Exception("Requested memory exceeds available memory")
        elif arg is "genome_url":
            data["genome"]["ncbi_url"] = getattr(args, arg)
        elif arg is "genome_file":
            if getattr(args, arg).split(".")[-1] in ["fasta", "fa", "fna"]:
                data["genome"]["file"] = getattr(args, arg)
            else:
                raise Exception("Provided reference file is not in an accepted format (.fasta, .fa, .fna).")

        with open(file, 'w') as yaml_file:
            yaml_file.write(yaml.dump(data, default_flow_style=False))
