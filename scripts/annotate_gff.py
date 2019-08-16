#!/usr/bin/env python

__authors__ = ["David Levy-Booth", "Parker Lloyd"]
__license__ = "GPL3"

# import argparse
#
# parser = argparse.ArgumentParser(description='gff annotation')
# parser.add_argument('gff_file', type=str, help='Input gff')
# parser.add_argument('blast_file', type=str, help='blast output (tab seperated)')
# parser.add_argument('new_gff', type=str, help='Output gff with annotation')
#
# args = parser.parse_args()
# gff_file=args.gff_file
# blast_file=args.blast_file
# new_gff=args.new_gff
gff_file=snakemake.input[0]
blast_file=snakemake.input[1]
new_gff=snakemake.output[0]

with open (blast_file, "r") as annotation:
    anotation_dict = {}
    for entry in annotation:
        entry = entry.split("\t")
        #print(entry)
        if entry: #test whether it is an empty line
            anotation_dict[entry[0]]=entry[1]
        else:
            continue

ofile = open(new_gff, "w")

with open (gff_file, "r") as mygff:
    for line in mygff:
        if not line.startswith ("#"):
            line_list = line.split(";")
            line_list = [i.split() for i in line_list]
            line_list = [item for sublist in line_list for item in sublist]
            gene = line_list[8].split("_")[1]
            contig = line_list[0]
            gene_id = contig+"_"+gene
            if line.rstrip().endswith(';'):
                new_line=line.rstrip()+"gene_id="+gene_id+";"
            else:
                new_line=line.rstrip()+";gene_id="+gene_id+";"
            if gene_id in anotation_dict:
                new_line=new_line+"product="+str(anotation_dict[gene_id])+";"
                ofile.write(new_line + "\n")
        else:
            ofile.write(line)

ofile.close()
