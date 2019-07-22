#!/usr/bin/env python

__authors__ = ["David Levy-Booth", "Parker Lloyd"]
__license__ = "GPL3"

# fasta_file=snakemake.input[0]
# blast_file=snakemake.input[1]
# new_fasta=snakemake.output[0]
fasta_file="/home/david/Liard/metaWRAP/Liard_Binning_refiner_outputs/BR_prodigal/Liard_36/Liard_36.clean.faa"
blast_file="/home/david/Liard/metaWRAP/Liard_Binning_refiner_outputs/BR_prodigal/Liard_36/Liard_36.dmnd.out.reduced"
new_fasta="/home/david/Liard/metaWRAP/Liard_Binning_refiner_outputs/BR_prodigal/Liard_36/Liard_36.clean.annotated.faa"

with open (blast_file, "r") as annotation:
    anotation_dict = {}
    for entry in annotation:
        entry = entry.split("\t")
        if entry: #test whether it is an empty line
            anotation_dict[entry[0]]=entry[1]
        else:
            continue

ofile = open(new_fasta, "w")

with open (fasta_file, "r") as myfasta:
    for line in myfasta:
        if line.startswith (">"):
            line = line[1:] # skip the ">" character
            line = line.split()
            if line[0] in anotation_dict:
                new_line = ">" + str(line[0]) + "\t" + str(anotation_dict[line[0]])
                ofile.write ( new_line + "\n")
            else:
                ofile.write ( ">" + line[0] + "\n")
        else:
            ofile.write(line)

ofile.close()
