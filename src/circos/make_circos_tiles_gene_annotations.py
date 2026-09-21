# make gene annotation input file for circos tiles

# File format:
#  chromosome_name start_pos end_pos
# color=color_names


# for use in circos tiles, extracted from the gff annotation
# TODO implement different colors for tiles in autosomes/X/Y

import os

input_gff_file = "/Users/miltr339/data_bioinfo/assemblies/chinensis_from_uppmax.fasta"
#input_gff_file = "/Users/milena/work/rc3.1_chinensis.gff"
outfile_name = "circos_gene_tiles.txt"

colorX = "color=acen"
colorY = "color=blue"
colorA = "color=lorange"

with open(input_gff_file, "r") as input_gff, open(outfile_name, "w") as outfile:
    for gff_line in input_gff.readlines()[1:]: #skip first line with gff3 header
        gff_line = gff_line.split("\t")
        if gff_line[2] == "gene": # only start and end coordinates of the complete gene, don't look at exons individually
            out_line = " ".join((gff_line[0], gff_line[3], gff_line[4])) #(()) because join takes tuples as input
            outfile.write(out_line+"\n")
        else:
            pass
        #break
print(os.getcwd()+"/"+outfile_name)
