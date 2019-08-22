import os
import yaml
import argparse
import psutil
import copy
from os import listdir
from os.path import isfile, join

# Define flags
parser = argparse.ArgumentParser(description='pipeline configuration')
parser.add_argument("-d", "--method", type = str, help = "Either featureCounts, htseq, or salmon.")
parser.add_argument("-a", "--aligner", type = str, help = "Either bt2, bbmap, or none.")
parser.add_argument("-q", "--trimmer", type = str, help = "Either trimmomatic, bbduk, or untrimmed.")
parser.add_argument("-s", "--assembler", type = str, help = "Either trinity or megahit.")
parser.add_argument("-t", "--threads", type = int, help = "The number of threads the pipeline is allowed to use.")
parser.add_argument("-m", "--max_memory", type = int, help = "The amount of memory provided to the assembler. Enter in byte format.")
parser.add_argument("-u", "--genome_url", type = str, help = "The NCBI url from which the reference genome can be downloaded.")
parser.add_argument("-f", "--genome_file", type = str, help = "The file containing the reference genome.")
parser.add_argument("-f", "--genes_faa", type = str, help = "Predicted genes in .faa format.")
parser.add_argument("-f", "--genes_fna", type = str, help = "Predicted genes in .fna format.")
parser.add_argument("-f", "--genes_gff", type = str, help = "Predicted genes in .gff format.")


args = parser.parse_args()

# Open and create an image of the Snakemake config
file = "config.yaml"
stream = open(file, "r")
original_data = yaml.load(stream, yaml.FullLoader)
new_data = copy.deepcopy(original_data)

# Iterate through all possible flags and update corresponding field in config image if flag is present and value is accepted
for arg in vars(args):
    if getattr(args, arg) is not None:
        if arg is "method":
            val = getattr(args, arg)
            if val in ["featureCounts", "htseq", "salmon"]:
                new_data[arg.upper()] = val
            else:
                raise Exception("Requested method is not an available option.")
        elif arg is "aligner":
            val = getattr(args, arg)
            if val in ["bt2", "bbmap", "none"]:
                new_data[arg.upper()] = val
            else:
                raise Exception("Requested aligner is not an available option.")
        elif arg is "trimmer":
            val = getattr(args, arg)
            if val in ["trimmomatic", "bbduk", "untrimmed"]:
                new_data[arg.upper()] = val
            else:
                raise Exception("Requested trimmer is not an available option.")
        elif arg is "assembler":
            val = getattr(args, arg)
            if val in ["trinity", "megahit"]:
                new_data[arg.upper()] = val
            else:
                raise Exception("Requested assembler is not an available option.")
        elif arg is "threads":
            val = getattr(args, arg)
            if val <= psutil.cpu_count():
                new_data[arg.upper()] = val
            else:
                raise Exception("Number of threads requested exceeds number of available cores.")
        elif arg is "max_memory":
            val = getattr(args, arg)
            if val <= psutil.virtual_memory()[1]:
                new_data[arg.upper()] = val
            else:
                raise Exception("Requested memory exceeds available memory")
        elif arg is "genome_url":
            new_data["genome"]["ncbi_url"] = getattr(args, arg)
        elif arg is "genome_file":
            if getattr(args, arg).split(".")[-1] in ["fasta", "fa", "fna"]:
                new_data["genome"]["ref_file"] = getattr(args, arg)
            else:
                raise Exception("Provided reference file is not in an accepted format (.fasta, .fa, .fna).")
        elif arg is "genes_faa":
            if getattr(args, arg).split(".")[-1] == "faa":
                new_data["genome"]["genes_faa_file"] = getattr(args, arg)
            else:
                raise Exception("Predicted genes file is not in an accepted format (.faa).")
        elif arg is "genes_fna":
            if getattr(args, arg).split(".")[-1] == "fna":
                new_data["genome"]["genes_fna_file"] = getattr(args, arg)
            else:
                raise Exception("Predicted genes file is not in an accepted format (.fna).")
        elif arg is "genes_gff":
            if getattr(args, arg).split(".")[-1] == "gff":
                new_data["genome"]["genes_gff_file"] = getattr(args, arg)
            else:
                raise Exception("Predicted genes file is not in an accepted format (.gff).")

# Rewrite config file if changes were made to image
if original_data != new_data:
    with open(file, 'w') as yaml_file:
        yaml_file.write(yaml.dump(new_data, default_flow_style=False))

# Close file stream
stream.close()
