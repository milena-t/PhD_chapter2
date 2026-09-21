### make circos karyotype file


## Formatting of the circos karyotype file:

# name, label, start and end position and a color
# chr - ID LABEL START END COLOR
# for example:
# chr - hs1 1 0 249250621 chr1
# chr - hs2 2 0 243199373 chr2
# chr - hs3 3 0 198022430 chr3

### based on the samtools faidx file, where columns 0 and 1 are the contig name and length




from Bio import SeqIO
import pandas as pd
import os
from collections import Counter
import sex_chromosomes

def fai_files(username="miltr339"):
    files_dir = f"/Users/{username}/work/chapter2/circos/karyotype/"
    outdict = {
        "A_obtectus" : f"{files_dir}A_obtectus.masked.fna.fai",
        "B_siliquastri" : f"{files_dir}B_siliquastri.masked.fna.fai",
        "B_varius" : f"{files_dir}B_varius.masked.fna.fai",
        "C_chinensis" : f"{files_dir}C_chinensis.masked.fna.fai",
        "C_maculatus" : f"{files_dir}C_maculatus.masked.fna.fai",
        "D_carinulata" : f"{files_dir}D_carinulata.masked.fna.fai",
        "D_sublineata" : f"{files_dir}D_sublineata.masked.fna.fai",
    }
    return outdict

def autosomes_lists():
    outdict = {
        "A_obtectus" : ["CAVLJG010000001.1","CAVLJG010000002.1","CAVLJG010000003.1","CAVLJG010000004.1","CAVLJG010000005.1","CAVLJG010000006.1","CAVLJG010000007.1","CAVLJG010000008.1","CAVLJG010000009.1"],
        "B_siliquastri" : ["1","2","3","4","5","6","7","8","9"],
        "B_varius" : ["OZ123443.1","OZ123444.1","OZ123445.1","OZ123446.1","OZ123447.1","OZ123448.1","OZ123449.1","OZ123450.1"],
        "C_chinensis" : [],
        "C_maculatus" : ["scaffold_1","scaffold_2","scaffold_3","scaffold_4","scaffold_5","scaffold_6","scaffold_7","scaffold_8","scaffold_9"],
        "D_carinulata" : ["NC_079461.1","NC_079462.1","NC_079463.1","NC_079464.1","NC_079465.1","NC_079466.1","NC_079467.1","NC_079468.1","NC_079469.1","NC_079470.1","NC_079471.1","NC_079472.1"],
        "D_sublineata" : ["NC_079474.1","NC_079475.1","NC_079476.1","NC_079477.1","NC_079478.1","NC_079479.1","NC_079480.1","NC_079481.1","NC_079482.1","NC_079483.1","NC_079484.1"],
    }
    return outdict

def make_karyotype_file(fai_filename,A_list, sex_chr_list,outfile_name, min_contig_len = 0):

    circos_cols= {
        "A" : "lorange",
        "X" : "acen",
        "Y" : "blue"}
    
    if A_list != []:
        min_contig_len = 0

    with open(fai_filename, "r") as fai_file, open(outfile_name, "w") as outfile:
        chr_num = 1
        for line in fai_file.readlines():
            contig,length,_,_,_ = line.strip().split("\t")

            if contig in A_list:
                col = circos_cols["A"]
                outfile_line = f"chr - {contig} {chr_num} 0 {length} {col}\n"
                outfile.write(outfile_line)
                chr_num += 1
                continue
            elif contig in sex_chr_list["X"]:
                col = circos_cols["X"]
                outfile_line = f"chr - {contig} {contig} 0 {length} {col}\n"
                outfile.write(outfile_line)
                continue
            elif contig in sex_chr_list["Y"]:
                col = circos_cols["Y"]
                outfile_line = f"chr - {contig} {contig} 0 {length} {col}\n"
                outfile.write(outfile_line)
                continue

            elif min_contig_len >0 and int(length)>min_contig_len:
                col = circos_cols["A"]
                outfile_line = f"chr - {contig} {contig} 0 {length} {col}\n"
                outfile.write(outfile_line)
                continue
            
    print(f"karyotype file written to: {outfile_name}")


if __name__=="__main__":
    username="miltr339"
    sex_chromosomes_dict = sex_chromosomes.get_contig_names()
    fai_dict = fai_files(username=username)
    autosomes_dict = autosomes_lists()
    data_dir = f"/Users/{username}/work/PhD_code/PhD_chapter2/data/circos"

    for species in sex_chromosomes_dict.keys():
        print(f"\n----------- {species} -----------")
        make_karyotype_file(fai_dict[species], autosomes_dict[species], sex_chr_list=sex_chromosomes_dict[species], outfile_name=f"{data_dir}/{species}_circos_karyotype.txt", min_contig_len = 7500000)