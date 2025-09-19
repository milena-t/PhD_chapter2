"""
MCScanX https://github.com/wyp1125/MCScanX requires a simplified annotation in bed format. This script takes a normal gff
annotation as input and returns a bed file suitable for MCScanX
the bed file is specificed like this (tab separated): 
chr#    starting_position       ending_position gene_ID
use like this: 
python3 make_bedfile_for_MCScanX.py annotation.gff
"""

import sys
import os

def modify_gff_for_MCScanX(gff_path:str):
    """ modify the file for MCScanX into a bedfile """
    outfile_bed_path = gff_path.split(".")[0]
    outfile_bed_path = f"{outfile_bed_path}_simplified.bed"

    with open(gff_path, "r") as gff_infile , open(outfile_bed_path, "w") as bed_outfile:
        for line in gff_infile.readlines():
            line = line.strip()
            if line[0] == "#":
                # exclude comment lines
                continue
            
            # parse GFF line
            contig,source,category_,start,stop,score,strandedness,frame,attributes=[c for c in line.split("\t") if len(c)>0]
            if category_ != "gene":
                continue
            gene_id = attributes.strip().split(";")[0].split("ID")[-1].replace("=","")
            
            bed_string = f"{contig}\t{start}\t{stop}\t{gene_id}\n"
            bed_outfile.write(bed_string)
    
    print(f"outfile written to {outfile_bed_path}")
            


if __name__ == "__main__":

    if len(sys.argv) == 1:
        annotation_file = "/Users/miltr339/work/native_annotations/a_obtectus_native_isoform_filtered.gff"
    elif len(sys.argv) == 2:
        annotation_file = sys.argv[2]
    else:
        raise RuntimeError(f"please specify one annotatin file as command line argument.\nyou have specified:\n{sys.argv}")

    if not os.access(annotation_file, os.R_OK):
        raise RuntimeError(f"{annotation_file} does not exist or is not readable")

    modify_gff_for_MCScanX(annotation_file)