# make the input file for circos that shows the blast hits

## output Format: (tab separated)
# source_contig_name	start_coordinate	end_coordinate	target_contig_name	start_coordinate	end_coordinate	color=color_name

# if some transcript IDs have "_1" at the end that match the annotations:
# for FILE in *species* ; do echo $FILE ; awk 'BEGIN{FS=OFS="\t"} {sub(/_1$/, "", $1)} 1' $FILE  > tmp && mv tmp $FILE ; done
# for FILE in *species* ; do echo $FILE ; awk 'BEGIN{FS=OFS="\t"} {sub(/_1$/, "", $2)} 1' $FILE  > tmp && mv tmp $FILE ; done

import pandas as pd
import os
import parse_gff as gff
import sex_chromosomes

def blast_paths_cmac_populations(username="miltr339"):
    filesdir = f"/Users/{username}/work/chapter2/circos/blast_outfiles"
    # /proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/protein_data/C_maculatus_populations
    # awk -F'\t' -v OFS='\t' '{ sub(/_1$/, "", $1); print }' C_maculatus_vs_C_maculatusC.blast > tmp.$$ && mv tmp.$$ C_maculatus_vs_C_maculatusC.blast
    # awk -F'\t' -v OFS='\t' '{ sub(/_1$/, "", $2); print }' C_maculatusC_vs_C_maculatus.blast > tmp.$$ && mv tmp.$$ C_maculatusC_vs_C_maculatus.blast

    outdir = {
        "China" : {
            "Lome" : f"{filesdir}/C_maculatusC_vs_C_maculatus.blast",
            "China" : f"{filesdir}/C_maculatusC_vs_C_maculatusC.blast",
        },
        "Lome" : {
            "China" : f"{filesdir}/C_maculatus_vs_C_maculatusC.blast",
        }
    }
    return outdir


def get_blast_paths(username="miltr339"):
    filedir = f"/Users/{username}/work/chapter2/circos/blast_outfiles/" 
    # /proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/protein_data/blastp_results
    outdir = {
        "A_obtectus" : {
            "A_obtectus" : f"{filedir}A_obtectus_vs_A_obtectus.blast",
            "B_siliquastri" : f"{filedir}A_obtectus_vs_B_siliquastri.blast",
            "B_varius" : f"{filedir}A_obtectus_vs_B_varius.blast",
            "C_chinensis" : f"{filedir}A_obtectus_vs_C_chinensis.blast",
            "C_maculatus" : f"{filedir}A_obtectus_vs_C_maculatus.blast",
            "D_carinulata" : f"{filedir}A_obtectus_vs_D_carinulata.blast",
            "D_sublineata" : f"{filedir}A_obtectus_vs_D_sublineata.blast",
            },
        "B_siliquastri" : {    
            "A_obtectus" : f"{filedir}B_siliquastri_vs_A_obtectus.blast",
            "B_siliquastri" : f"{filedir}B_siliquastri_vs_B_siliquastri.blast",
            "B_varius" : f"{filedir}B_siliquastri_vs_B_varius.blast",
            "C_chinensis" : f"{filedir}B_siliquastri_vs_C_chinensis.blast",
            "C_maculatus" : f"{filedir}B_siliquastri_vs_C_maculatus.blast",
            "D_carinulata" : f"{filedir}B_siliquastri_vs_D_carinulata.blast",
            "D_sublineata" : f"{filedir}B_siliquastri_vs_D_sublineata.blast",
            },
        "B_varius" : {    
            "A_obtectus" : f"{filedir}B_varius_vs_A_obtectus.blast",
            "B_siliquastri" : f"{filedir}B_varius_vs_B_siliquastri.blast",
            "B_varius" : f"{filedir}B_varius_vs_B_varius.blast",
            "C_chinensis" : f"{filedir}B_varius_vs_C_chinensis.blast",
            "C_maculatus" : f"{filedir}B_varius_vs_C_maculatus.blast",
            "D_carinulata" : f"{filedir}B_varius_vs_D_carinulata.blast",
            "D_sublineata" : f"{filedir}B_varius_vs_D_sublineata.blast",
            },
        "C_chinensis" : {    
            "A_obtectus" : f"{filedir}C_chinensis_vs_A_obtectus.blast",
            "B_siliquastri" : f"{filedir}C_chinensis_vs_B_siliquastri.blast",
            "B_varius" : f"{filedir}C_chinensis_vs_B_varius.blast",
            "C_chinensis" : f"{filedir}C_chinensis_vs_C_chinensis.blast",
            "C_maculatus" : f"{filedir}C_chinensis_vs_C_maculatus.blast",
            "D_carinulata" : f"{filedir}C_chinensis_vs_D_carinulata.blast",
            "D_sublineata" : f"{filedir}C_chinensis_vs_D_sublineata.blast",
            },
        "C_maculatus" : {    
            "A_obtectus" : f"{filedir}C_maculatus_vs_A_obtectus.blast",
            "B_siliquastri" : f"{filedir}C_maculatus_vs_B_siliquastri.blast",
            "B_varius" : f"{filedir}C_maculatus_vs_B_varius.blast",
            "C_chinensis" : f"{filedir}C_maculatus_vs_C_chinensis.blast",
            "C_maculatus" : f"{filedir}C_maculatus_vs_C_maculatus.blast",
            "D_carinulata" : f"{filedir}C_maculatus_vs_D_carinulata.blast",
            "D_sublineata" : f"{filedir}C_maculatus_vs_D_sublineata.blast",
            },
        "D_carinulata" : {    
            "A_obtectus" : f"{filedir}D_carinulata_vs_A_obtectus.blast",
            "B_siliquastri" : f"{filedir}D_carinulata_vs_B_siliquastri.blast",
            "B_varius" : f"{filedir}D_carinulata_vs_B_varius.blast",
            "C_chinensis" : f"{filedir}D_carinulata_vs_C_chinensis.blast",
            "C_maculatus" : f"{filedir}D_carinulata_vs_C_maculatus.blast",
            "D_carinulata" : f"{filedir}D_carinulata_vs_D_carinulata.blast",
            "D_sublineata" : f"{filedir}D_carinulata_vs_D_sublineata.blast",
            },
        "D_sublineata" : {    
            "A_obtectus" : f"{filedir}D_sublineata_vs_A_obtectus.blast",
            "B_siliquastri" : f"{filedir}D_sublineata_vs_B_siliquastri.blast",
            "B_varius" : f"{filedir}D_sublineata_vs_B_varius.blast",
            "C_chinensis" : f"{filedir}D_sublineata_vs_C_chinensis.blast",
            "C_maculatus" : f"{filedir}D_sublineata_vs_C_maculatus.blast",
            "D_carinulata" : f"{filedir}D_sublineata_vs_D_carinulata.blast",
            "D_sublineata" : f"{filedir}D_sublineata_vs_D_sublineata.blast",
            },
    }
    return outdir



def annotation_paths_cmac_populations(username="milena"):
    filesdir = f"/Users/{username}/work/chapter2/native_annotations"
    outdict = {
        "China" : f"{filesdir}/C_maculatusC.gff",
        "Lome" : f"{filesdir}/C_maculatus.gff",
    }
    return outdict

def get_annotation_paths(username="miltr339"):
    dirpath = f"/Users/{username}/work/chapter2/native_annotations/"
    outdict = {
        "A_obtectus" : f"{dirpath}A_obtectus.gff",
        "B_siliquastri" : f"{dirpath}B_siliquastri.gff",
        "B_varius" : f"{dirpath}B_varius.gff",
        "C_chinensis" : f"{dirpath}C_chinensis.gff",
        "C_maculatus" : f"{dirpath}C_maculatus.gff",
        "D_carinulata" : f"{dirpath}D_carinulata.gff",
        "D_sublineata" : f"{dirpath}D_sublineata.gff",
    }
    return outdict


def make_circos_hits_file(annotation_file1, annotation_file2, blast_outfile, sex_chr_contigs_dict1,sex_chr_contigs_dict2, circos_outfile_name, min_seq_ident=90, max_seq_ident = 200):

    try:
        gff1_dict = gff.parse_gff3_general(annotation_file1)
    except:
        gff1_dict = gff.parse_gff3_general(annotation_file1, gtf=True)
    try:
        gff2_dict = gff.parse_gff3_general(annotation_file2)
    except:
        gff2_dict = gff.parse_gff3_general(annotation_file2, gtf=True)

    colors = {
        "X" : "acen",
        "Y" : "blue",
        "A" : "lorange"
    }
    # blast_outfmt6_headers = ["qseqid", "rseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    circos_outfile_name_X = circos_outfile_name.replace(".txt", "_X.txt")
    circos_outfile_name_Y = circos_outfile_name.replace(".txt", "_Y.txt")
    with open(blast_outfile, "r") as filtered_outfile, open(circos_outfile_name, "w") as circos_outfile, open(circos_outfile_name_X, "w") as circos_outfile_X, open(circos_outfile_name_Y, "w") as circos_outfile_Y:
        for line_str in filtered_outfile:
            
            transcriptID1,transcriptID2,seq_ident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore = line_str.strip().split("\t")
            if float(seq_ident) < min_seq_ident:
                continue
            if float(seq_ident) >= max_seq_ident:
                continue
            
            try:
                transcript1 = gff1_dict[transcriptID1]
            except:
                raise RuntimeError(f"{transcriptID1} not found in {annotation_file1}!")
            try:
                transcript2 = gff2_dict[transcriptID2]
            except:
                raise RuntimeError(f"{transcriptID2} not found in {annotation_file2}!")

            contig1 = transcript1.contig
            contig2 = transcript2.contig

            if contig1 in sex_chr_contigs_dict1["X"]:
                color = colors["X"]
            elif contig1 in sex_chr_contigs_dict1["Y"]:
                color = colors["Y"]
            else:
                color = colors["A"]
            
            start1=transcript1.start
            end1=transcript1.end
            start2=transcript2.start
            end2=transcript2.end
            
            circos = f"{contig1} {start1} {end1} {contig2} {start2} {end2} color={color}"
            circos_outfile.write(circos+"\n")

            for chr in ["X","Y"]:
                other_chr = "Y"
                if chr == "Y":
                    other_chr = "X"
                if contig1 in sex_chr_contigs_dict1[chr] and contig2 in sex_chr_contigs_dict2[chr]:
                    color=colors[chr]
                elif (contig1 in sex_chr_contigs_dict1[chr] and contig2 in sex_chr_contigs_dict2[other_chr]) or (contig1 in sex_chr_contigs_dict1[other_chr] and contig2 in sex_chr_contigs_dict2[chr]):
                    color=colors[other_chr]
                elif (contig1 in sex_chr_contigs_dict1[chr] and contig2 not in sex_chr_contigs_dict2[chr]) or (contig1 not in sex_chr_contigs_dict1[chr] and contig2 in sex_chr_contigs_dict2[chr]):
                    color=colors["A"]
                else:
                    continue

                circos = f"{contig1} {start1} {end1} {contig2} {start2} {end2} color={color}"
                if chr=="X":
                    circos_outfile_X.write(circos+"\n")
                elif chr=="Y":
                    circos_outfile_Y.write(circos+"\n")

            #break
    print(f"The circos outfile is here: \n  - {circos_outfile_name}\n  - {circos_outfile_name_X}\n  - {circos_outfile_name_Y}")


if __name__ == "__main__":

    username="miltr339"

    blast_outfiles_dict = get_blast_paths(username=username)
    annotations_dict = get_annotation_paths(username=username)
    sex_chromosomes_dict = sex_chromosomes.get_contig_names()
    species_list = list(blast_outfiles_dict.keys())
    data_dir = f"/Users/{username}/work/PhD_code/PhD_chapter2/data/circos/"

    if True:
        for species1 in species_list:
            for species2 in species_list:
                if species1==species2:
                    max_seq_ident=100 # exclude self-hits for self-blast
                else:
                    max_seq_ident=200

                print(f"\n ------------ {species1} vs. {species2} ------------")
                make_circos_hits_file(annotation_file1=annotations_dict[species1],
                    annotation_file2=annotations_dict[species2],
                    blast_outfile=blast_outfiles_dict[species1][species2],
                    sex_chr_contigs_dict1=sex_chromosomes_dict[species1],
                    sex_chr_contigs_dict2=sex_chromosomes_dict[species2],
                    min_seq_ident=90,max_seq_ident=max_seq_ident,
                    circos_outfile_name=f"{data_dir}circos_links_{species1}_vs_{species2}.txt")

    if False:
        blast_outfiles_dict = blast_paths_cmac_populations(username=username)
        annotations_dict = annotation_paths_cmac_populations(username=username)
        sex_chromosomes_dict = sex_chromosomes.Cmac_S_L_nonscaffolded_contig_names()
        species_list = list(blast_outfiles_dict.keys())
        data_dir = f"/Users/{username}/work/PhD_code/PhD_chapter2/data/circos/Cmac_populations/"

        for species1 in species_list:
            for species2 in species_list:
                if species1=="Lome" and species2=="Lome" :
                    continue

                if species1==species2:
                    max_seq_ident=100 # exclude self-hits for self-blast
                else:
                    max_seq_ident=200

                # if species1=="China" and species2=="China" :
                #     print(f"------------------ {sex_chromosomes_dict[species1]}")
                #     print(f"------------------ {sex_chromosomes_dict[species2]}")

                print(f"\n ------------ {species1} vs. {species2} ------------")
                make_circos_hits_file(annotation_file1=annotations_dict[species1],
                    annotation_file2=annotations_dict[species2],
                    blast_outfile=blast_outfiles_dict[species1][species2],
                    sex_chr_contigs_dict1=sex_chromosomes_dict[species1],
                    sex_chr_contigs_dict2=sex_chromosomes_dict[species2],
                    min_seq_ident=90,max_seq_ident=max_seq_ident,
                    circos_outfile_name=f"{data_dir}circos_links_{species1}_vs_{species2}.txt")