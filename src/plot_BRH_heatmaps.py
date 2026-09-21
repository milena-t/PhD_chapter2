"""
Plot species-pair heatmaps of BRH chromosome locations to see if there are A/X/Y translocations
"""

import numpy as np
import scipy.stats as sts
import math
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

def species_pair_paths(username="miltr339"):
    dirpath = f"/Users/{username}/work/PhD_code/PhD_chapter2/data/BRH_orthologs/"
    out_dict = {
        "A_obtectus" : {
            "B_siliquastri" : f"{dirpath}A_obtectus_vs_B_siliquastri_BRHs.tsv",
            "B_varius" : f"{dirpath}A_obtectus_vs_B_varius_BRHs.tsv",
            "C_chinensis" : f"{dirpath}A_obtectus_vs_C_chinensis_BRHs.tsv",
            "C_maculatus" : f"{dirpath}A_obtectus_vs_C_maculatus_BRHs.tsv",
            "D_carinulata" : f"{dirpath}A_obtectus_vs_D_carinulata_BRHs.tsv",
            "D_sublineata" : f"{dirpath}A_obtectus_vs_D_sublineata_BRHs.tsv",
            },
        "B_siliquastri" : {
            "A_obtectus" : f"{dirpath}B_siliquastri_vs_A_obtectus_BRHs.tsv",
            "B_varius" : f"{dirpath}B_siliquastri_vs_B_varius_BRHs.tsv",
            "C_chinensis" : f"{dirpath}B_siliquastri_vs_C_chinensis_BRHs.tsv",
            "C_maculatus" : f"{dirpath}B_siliquastri_vs_C_maculatus_BRHs.tsv",
            "D_carinulata" : f"{dirpath}B_siliquastri_vs_D_carinulata_BRHs.tsv",
            "D_sublineata" : f"{dirpath}B_siliquastri_vs_D_sublineata_BRHs.tsv",
            },
        "B_varius" : {
            "A_obtectus" : f"{dirpath}B_varius_vs_A_obtectus_BRHs.tsv",
            "B_siliquastri" : f"{dirpath}B_varius_vs_B_siliquastri_BRHs.tsv",
            "C_chinensis" : f"{dirpath}B_varius_vs_C_chinensis_BRHs.tsv",
            "C_maculatus" : f"{dirpath}B_varius_vs_C_maculatus_BRHs.tsv",
            "D_carinulata" : f"{dirpath}B_varius_vs_D_carinulata_BRHs.tsv",
            "D_sublineata" : f"{dirpath}B_varius_vs_D_sublineata_BRHs.tsv",
            },
        "C_chinensis" : {
            "A_obtectus" : f"{dirpath}C_chinensis_vs_A_obtectus_BRHs.tsv",
            "B_siliquastri" : f"{dirpath}C_chinensis_vs_B_siliquastri_BRHs.tsv",
            "B_varius" : f"{dirpath}C_chinensis_vs_B_varius_BRHs.tsv",
            "C_maculatus" : f"{dirpath}C_chinensis_vs_C_maculatus_BRHs.tsv",
            "D_carinulata" : f"{dirpath}C_chinensis_vs_D_carinulata_BRHs.tsv",
            "D_sublineata" : f"{dirpath}C_chinensis_vs_D_sublineata_BRHs.tsv",
            },
        "C_maculatus" : {
            "A_obtectus" : f"{dirpath}C_maculatus_vs_A_obtectus_BRHs.tsv",
            "B_siliquastri" : f"{dirpath}C_maculatus_vs_B_siliquastri_BRHs.tsv",
            "B_varius" : f"{dirpath}C_maculatus_vs_B_varius_BRHs.tsv",
            "C_chinensis" : f"{dirpath}C_maculatus_vs_C_chinensis_BRHs.tsv",
            "D_carinulata" : f"{dirpath}C_maculatus_vs_D_carinulata_BRHs.tsv",
            "D_sublineata" : f"{dirpath}C_maculatus_vs_D_sublineata_BRHs.tsv",
            },
        "D_carinulata" : {
            "A_obtectus" : f"{dirpath}D_carinulata_vs_A_obtectus_BRHs.tsv",
            "B_siliquastri" : f"{dirpath}D_carinulata_vs_B_siliquastri_BRHs.tsv",
            "B_varius" : f"{dirpath}D_carinulata_vs_B_varius_BRHs.tsv",
            "C_chinensis" : f"{dirpath}D_carinulata_vs_C_chinensis_BRHs.tsv",
            "C_maculatus" : f"{dirpath}D_carinulata_vs_C_maculatus_BRHs.tsv",
            "D_sublineata" : f"{dirpath}D_carinulata_vs_D_sublineata_BRHs.tsv",
            },
        "D_sublineata" : {
            "A_obtectus" : f"{dirpath}D_sublineata_vs_A_obtectus_BRHs.tsv",
            "B_siliquastri" : f"{dirpath}D_sublineata_vs_B_siliquastri_BRHs.tsv",
            "B_varius" : f"{dirpath}D_sublineata_vs_B_varius_BRHs.tsv",
            "C_chinensis" : f"{dirpath}D_sublineata_vs_C_chinensis_BRHs.tsv",
            "C_maculatus" : f"{dirpath}D_sublineata_vs_C_maculatus_BRHs.tsv",
            "D_carinulata" : f"{dirpath}D_sublineata_vs_D_carinulata_BRHs.tsv",
            },
    }
    return out_dict

def make_BRH_table_chr_counts(table_path, chr_list = ["A", "X", "Y"]):
    count_array = np.zeros((len(chr_list),len(chr_list)))
    df_full = pd.read_csv(table_path, sep="\t")
    for i, chr1 in enumerate(chr_list):
        df_filt1 = df_full[df_full["chromosome"]==chr1]
        for j, chr2 in enumerate(chr_list):
            df_filt2 = df_filt1[df_filt1["chromosome.1"]==chr2]
            # print(f"{chr1} + {chr2} = {len(df_filt2)}")
            count_array[i][j] = len(df_filt2)
    return count_array.astype(int)
    

def heatmap_show_counts(row,col,chr_list,count_array,axes,fs):
    for i in range(len(chr_list)):
        for j in range(len(chr_list)):
            try:
                count = count_array[i, j]
                text = axes[row,col].text(j, i, f"{count}",ha="center", va="center", color="w", fontsize = fs*0.75)
            except:
                continue



def plot_pairwise_heatmaps(pair_tables:dict, plot_filename:str="heatmap.png", chr_list=["A", "X", "Y"], log_colors = True, species_order = []):

    plt.rcParams['text.usetex'] = True # use \\textit{{{}}} for species names
    plt.rcParams['text.latex.preamble'] = r'\usepackage{sfmath} \renewcommand{\familydefault}{\sfdefault}'
    plt.rcParams['font.family'] = 'sans-serif'

    if species_order==[]:
        species_list = list(pair_tables.keys())
    else:
        species_list = species_order
    species_count = len(species_list)
    cols = species_count
    rows = cols
    if rows>2:
        fig, axes = plt.subplots(rows, cols, figsize=(25, 25)) # for more than three rows
    else:
        fig, axes = plt.subplots(rows, cols, figsize=(15, 10)) # for more than three rows
    
    fs = 25

    for row, species1 in enumerate(species_list):
        for col, species2 in enumerate(species_list):
            if species1==species2:
                # plot species name
                species1_=species1.replace("_",". ")
                axes[row,col].axis('off')
                axes[row,col].text(0.15,0.4,f"\\textit{{{species1_}}}", fontsize = fs*1.25)# rotation = 90, 
                continue

            count_array = make_BRH_table_chr_counts(pair_tables[species1][species2], chr_list=chr_list)
            # print(f"\n {species1} + {species2}:\n{count_array}")
            print(f"\n {species1} vs. {species2}")

            cmap = LinearSegmentedColormap.from_list("blue_to_green",["#38ABD2", "#008F69"])
            # cbarlabel="number of pairwise 1-to-1 ortholgs"
            if log_colors:
                count_array_ = np.zeros((len(chr_list),len(chr_list)))
                for i in range(len(chr_list)):
                    for j in range(len(chr_list)):
                        count_array_[i, j] =  math.log(count_array[i, j]+1)
                im = axes[row,col].imshow(count_array_, cmap=cmap)
            else:
                im = axes[row,col].imshow(count_array, cmap=cmap)

            heatmap_show_counts(row=row,col=col,chr_list=chr_list,count_array=count_array,axes=axes,fs=fs)

            axes[row,col].set_xticks(range(len(chr_list)), labels=chr_list, fontsize = fs) # ha="right", rotation_mode="anchor", rotation=90
            axes[row,col].set_yticks(range(len(chr_list)), labels=chr_list, fontsize = fs)
            
            ## the array coordinates are [reference, query]
            # axes[row,col].set_xlabel("query", fontsize=fs)
            # axes[row,col].set_ylabel("reference", fontsize=fs)

    # plt.title(label="1-to-1 ortholog chromosome locations", fontsize = fs*1.3)
    
    plt.tight_layout()
    # plt.show()
    plt.savefig(plot_filename, dpi = 300, transparent = True)
    print(f"figure saved here: {plot_filename}")



if __name__ == "__main__":
    username = "miltr339"
    pair_tables = species_pair_paths(username=username)

    heatmap_name=f"/Users/{username}/work/PhD_code/PhD_chapter2/data/BRH_orthologs/counts_heatmap.png"
    species_order = [
"C_maculatus",
"C_chinensis",
"A_obtectus",
"B_siliquastri",
"B_varius",
"D_carinulata",
"D_sublineata"]
    plot_pairwise_heatmaps(pair_tables=pair_tables, plot_filename=heatmap_name, species_order=species_order)