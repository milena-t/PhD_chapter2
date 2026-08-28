"""
Plot a heatmap for all pairwise genome alignments. I will be inspired by this: https://www.nature.com/articles/s41586-024-07473-2/figures/1
and they do: "percentage of reference, as shown on the x axis, covered by the query, as shown on the y axis"
"""

import sex_chromosomes

import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

def fasta_indices(username="milena"):
    """
    read the assembly index files and make dicts of {
        species : {contig : length, contig : length , ...}
    }
    rsync -azP "milenatr@pelle.uppmax.uu.se:/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/assemblies/*fai" /Users/miltr339/work/chapter2/wga/
    """
    dirname = f"/Users/{username}/work/chapter2/wga/"
    filenames_dict = {
        "A_obtectus" : f"{dirname}A_obtectus.masked.fna.fai",
        "B_siliquastri" : f"{dirname}B_siliquastri.masked.fna.fai",
        "B_varius" : f"{dirname}B_varius.masked.fna.fai",
        "C_chinensis" : f"{dirname}C_chinensis.masked.fna.fai",
        "C_maculatus" : f"{dirname}C_maculatus.masked.fna.fai",
        "D_carinulata" : f"{dirname}D_carinulata.masked.fna.fai",
        "D_sublineata" : f"{dirname}D_sublineata.masked.fna.fai",
    }
    def read_faidx(infile_path):
        df = pd.read_csv(infile_path, sep="\t", header=None)
        return {contig : length for contig,length in zip(df.iloc[:,0],df.iloc[:,1])}

    out_dict = {species : read_faidx(filename) for species,filename in filenames_dict.items()}
    
    dirname=f"/Users/{username}/work/chapter2/wga/Cmac_populations/"
    filenames_dict = {
        "China" : f"{dirname}Cmac_china.fasta.fai",
        "Lome_Yl" : f"{dirname}Lome_Yl.fasta.fai",
        "Lome_Ys" : f"{dirname}Lome_Ys.fasta.fai",
    }
    out_dict_cmac = {species : read_faidx(filename) for species,filename in filenames_dict.items()}

    return out_dict,out_dict_cmac

def get_aln_coord_files(username="milena"):
    dirname = f"/Users/{username}/work/chapter2/wga/"
    out_dict = {
        "A_obtectus" : {
            "B_siliquastri" : f"{dirname}ref_A_obtectus.masked_query_B_siliquastri.masked.delta_filtered_coords",
            "B_varius" : f"{dirname}ref_A_obtectus.masked_query_B_varius.masked.delta_filtered_coords",
            "C_chinensis" : f"{dirname}ref_A_obtectus.masked_query_C_chinensis.masked.delta_filtered_coords",
            "C_maculatus" : f"{dirname}ref_A_obtectus.masked_query_C_maculatus.masked.delta_filtered_coords",
            "D_carinulata" : f"{dirname}ref_A_obtectus.masked_query_D_carinulata.masked.delta_filtered_coords",
            "D_sublineata" : f"{dirname}ref_A_obtectus.masked_query_D_sublineata.masked.delta_filtered_coords",
        },
        "B_siliquastri" : {
            "A_obtectus" : f"{dirname}ref_B_siliquastri.masked_query_A_obtectus.masked.delta_filtered_coords",
            "B_varius" : f"{dirname}ref_B_siliquastri.masked_query_B_varius.masked.delta_filtered_coords",
            "C_chinensis" : f"{dirname}ref_B_siliquastri.masked_query_C_chinensis.masked.delta_filtered_coords",
            "C_maculatus" : f"{dirname}ref_B_siliquastri.masked_query_C_maculatus.masked.delta_filtered_coords",
            "D_carinulata" : f"{dirname}ref_B_siliquastri.masked_query_D_carinulata.masked.delta_filtered_coords",
            "D_sublineata" : f"{dirname}ref_B_siliquastri.masked_query_D_sublineata.masked.delta_filtered_coords",
        },
        "B_varius" : {
            "A_obtectus" : f"{dirname}ref_B_varius.masked_query_A_obtectus.masked.delta_filtered_coords",
            "B_siliquastri" : f"{dirname}ref_B_varius.masked_query_B_siliquastri.masked.delta_filtered_coords",
            "C_chinensis" : f"{dirname}ref_B_varius.masked_query_C_chinensis.masked.delta_filtered_coords",
            "C_maculatus" : f"{dirname}ref_B_varius.masked_query_C_maculatus.masked.delta_filtered_coords",
            "D_carinulata" : f"{dirname}ref_B_varius.masked_query_D_carinulata.masked.delta_filtered_coords",
            "D_sublineata" : f"{dirname}ref_B_varius.masked_query_D_sublineata.masked.delta_filtered_coords",

        },
        "C_chinensis" : {
            "A_obtectus" : f"{dirname}ref_C_chinensis.masked_query_A_obtectus.masked.delta_filtered_coords",
            "B_siliquastri" : f"{dirname}ref_C_chinensis.masked_query_B_siliquastri.masked.delta_filtered_coords",
            "B_varius" : f"{dirname}ref_C_chinensis.masked_query_B_varius.masked.delta_filtered_coords",
            "C_maculatus" : f"{dirname}ref_C_chinensis.masked_query_C_maculatus.masked.delta_filtered_coords",
            "D_carinulata" : f"{dirname}ref_C_chinensis.masked_query_D_carinulata.masked.delta_filtered_coords",
            "D_sublineata" : f"{dirname}ref_C_chinensis.masked_query_D_sublineata.masked.delta_filtered_coords",
        },
        "C_maculatus" : {
            "A_obtectus" : f"{dirname}ref_C_maculatus.masked_query_A_obtectus.masked.delta_filtered_coords",
            "B_siliquastri" : f"{dirname}ref_C_maculatus.masked_query_B_siliquastri.masked.delta_filtered_coords",
            "B_varius" : f"{dirname}ref_C_maculatus.masked_query_B_varius.masked.delta_filtered_coords",
            "C_chinensis" : f"{dirname}ref_C_maculatus.masked_query_C_chinensis.masked.delta_filtered_coords",
            "D_carinulata" : f"{dirname}ref_C_maculatus.masked_query_D_carinulata.masked.delta_filtered_coords",
            "D_sublineata" : f"{dirname}ref_C_maculatus.masked_query_D_sublineata.masked.delta_filtered_coords",
        },
        "D_carinulata" : {
            "A_obtectus" : f"{dirname}ref_D_carinulata.masked_query_A_obtectus.masked.delta_filtered_coords",
            "B_siliquastri" : f"{dirname}ref_D_carinulata.masked_query_B_siliquastri.masked.delta_filtered_coords",
            "B_varius" : f"{dirname}ref_D_carinulata.masked_query_B_varius.masked.delta_filtered_coords",
            "C_chinensis" : f"{dirname}ref_D_carinulata.masked_query_C_chinensis.masked.delta_filtered_coords",
            "C_maculatus" : f"{dirname}ref_D_carinulata.masked_query_C_maculatus.masked.delta_filtered_coords",
            "D_sublineata" : f"{dirname}ref_D_carinulata.masked_query_D_sublineata.masked.delta_filtered_coords",
        },
        "D_sublineata" : {
            "A_obtectus" : f"{dirname}ref_D_sublineata.masked_query_A_obtectus.masked.delta_filtered_coords",
            "B_siliquastri" : f"{dirname}ref_D_sublineata.masked_query_B_siliquastri.masked.delta_filtered_coords",
            "B_varius" : f"{dirname}ref_D_sublineata.masked_query_B_varius.masked.delta_filtered_coords",
            "C_chinensis" : f"{dirname}ref_D_sublineata.masked_query_C_chinensis.masked.delta_filtered_coords",
            "C_maculatus" : f"{dirname}ref_D_sublineata.masked_query_C_maculatus.masked.delta_filtered_coords",
            "D_carinulata" : f"{dirname}ref_D_sublineata.masked_query_D_carinulata.masked.delta_filtered_coords",
        }
    }

    dirname=f"/Users/{username}/work/chapter2/wga/Cmac_populations/"
    out_dict_cmac = {
        "China" : {
            "Lome_Yl" : f"{dirname}ref_Cmac_china_query_Lome_Yl.delta_filtered_coords",
            "Lome_Ys" : f"{dirname}ref_Cmac_china_query_Lome_Ys.delta_filtered_coords",
        },
        "Lome_Yl" : {
            "China" : f"{dirname}ref_Lome_Yl_query_Cmac_china.delta_filtered_coords",
            "Lome_Ys" : f"{dirname}ref_Lome_Yl_query_Lome_Ys.delta_filtered_coords",
        },
        "Lome_Ys" : {
            "China" : f"{dirname}ref_Lome_Ys_query_Cmac_china.delta_filtered_coords",
            "Lome_Yl" : f"{dirname}ref_Lome_Ys_query_Lome_Yl.delta_filtered_coords",
        },
    }

    return out_dict,out_dict_cmac

def count_ref_aln_coordinates(aln_filtered_coords_path, ref_faidx_dict, contig_list = [], verbose=False):
    """
    get the percentage of bases of each reference contig that are covered by the query. 
    full contig length is determined from the faidx file
    The first two column in the alignment coordinates are start and end of the alignment in the reference contig
    if given a contig list, the percentage is the coverage across all contigs
    """
    out_cov_bases = {}
    with open(aln_filtered_coords_path, "r") as aln_filtered_coords:
        for line in aln_filtered_coords.readlines():
            refstart,refend,qstart,qend,reflen,qlen,percident,refcontig,qcontig = line.strip().split("\t")                
            if refcontig in out_cov_bases:
                out_cov_bases[refcontig] += int(reflen)
            else:
                out_cov_bases[refcontig] = int(reflen)

    if contig_list == []:
        out_perc_bases = {contig : num_bases / ref_faidx_dict[contig] for contig,num_bases in out_cov_bases.items()}
        return out_perc_bases
    else:
        all_len = 0
        cov_len = 0
        miss_cont = []
        for contig in contig_list:
            if contig in out_cov_bases:
                cov_len += out_cov_bases[contig]
                all_len += ref_faidx_dict[contig]
            else:
                miss_cont.append(contig)
        if len(miss_cont) > 0 and verbose:
            uncov_len = sum([ref_faidx_dict[contig] for contig in miss_cont])
            print(f"\t{len(miss_cont)} contigs ({uncov_len} bp) out of {len(contig_list)} from input list not in alignment: {aln_filtered_coords_path} \n\t{miss_cont}")
        if all_len > 0:
            return cov_len / all_len
        else:
            return 0


def make_array_for_heatmap(aln_coord_files, sex_chr_dict, faidx_dicts, chr = "X",verbose=False, log_transform=False):
    """
    make np.array to plot heatmap
    """
    species_list = list(aln_coord_files.keys())
    species_count =len(species_list)
    species_index = {species : i for i,species in enumerate(species_list)}
    perc_overlap = np.empty((species_count,species_count))
    perc_overlap[:] = np.nan

    for ref_species,query_dict in aln_coord_files.items():
        if verbose:
            print(f"==================== {ref_species} ====================")
        index_1 = species_index[ref_species]
        ref_X = sex_chr_dict[ref_species][chr]
        for query_species,aln_coords_path in query_dict.items():
            index_2 = species_index[query_species]
            quot_bases = count_ref_aln_coordinates(aln_filtered_coords_path=aln_coords_path, ref_faidx_dict=faidx_dicts[ref_species], contig_list=ref_X)
            perc_bases = quot_bases*100
            if log_transform and perc_bases!= 0:
                perc_bases = math.log(perc_bases)
                
            perc_overlap[index_1, index_2] = perc_bases
            if verbose and not log_transform:
                print(f"{query_species} : {perc_bases:.4f} %")
    
    return perc_overlap,species_list



def plot_heatmap(counts_array, species_list, filename = "BRH_orthologs_heatmap.png", title = f"pairwise orthologs counts", log_data = False):
    """
    plot the heatmap created in make_array_for_heatmap()
    """

    fs = 40
    aspect_ratio = 18 / 18
    height_pixels = 2000  # Height in pixels
    width_pixels = int(height_pixels * aspect_ratio)  # Width in pixels
    fig = plt.figure(figsize=(width_pixels / 100, height_pixels / 100), dpi=100)
    ax = fig.add_subplot(111)

    # cmap=mpl.colormaps["rainbow"]
    cmap = LinearSegmentedColormap.from_list("red_to_orange",["#b82946", "#F2933A"])
    # cbarlabel="number of pairwise 1-to-1 ortholgs"
    im = ax.imshow(counts_array, cmap=cmap)

    # create text annotations
    for i in range(len(species_list)):
        for j in range(len(species_list)):
            if i==j:
                continue
            try:
                count = counts_array[i, j]
                if log_data: # back-transform to actual percentage, only log by color
                    count = math.exp(count)
                text = ax.text(j, i, f"{count:.1f}%",ha="center", va="center", color="w", fontsize = fs*0.75)
            except:
                continue

    species_list = [species.replace("_", ". ") for species in species_list]
    # Show all ticks and label them with the respective list entries
    ax.set_xticks(range(len(species_list)), labels=species_list,rotation=45, ha="right", rotation_mode="anchor", fontsize = fs)
    ax.set_yticks(range(len(species_list)), labels=species_list, fontsize = fs)
    plt.title(label=title, fontsize = fs*1.3)
    
    plt.tight_layout()
    # plt.show()
    plt.savefig(filename, dpi = 300, transparent = True)
    print(f"figure saved here: {filename}")


if __name__ == "__main__":
    username = "miltr339"
    faidx_dicts,faidx_dicts_cmac = fasta_indices(username=username)
    aln_coord_files,aln_coord_files_cmac = get_aln_coord_files(username=username)
    outdir_wga = f"/Users/{username}/work/PhD_code/PhD_chapter2/data/pairwise_wga"

    # when percentages are really low, log-transform the heatmap colors to still see the variation
    logtr = True
    if logtr:
        log_text = "log"
    else:
        log_text = ""

    if False:
        sex_chr_dict = sex_chromosomes.get_contig_names()
        for chr in ["X","Y"]:
            perc_X_overlap,species_list = make_array_for_heatmap(aln_coord_files=aln_coord_files, sex_chr_dict=sex_chr_dict, faidx_dicts=faidx_dicts, chr=chr,verbose=False, log_transform=logtr)
            print(perc_X_overlap)
            plot_heatmap(counts_array=perc_X_overlap, species_list=species_list, log_data=logtr,
            filename = f"{outdir_wga}/{chr}_chr_{log_text}_alignment_coverage_heatmap.png", title = f"{chr}-Chromosome aln. coverage")

    if True:
        sex_chr_dict = sex_chromosomes.Cmac_S_L_nonscaffolded_contig_names()
        for chr in ["X","Y"]:
            perc_X_overlap,species_list = make_array_for_heatmap(aln_coord_files=aln_coord_files_cmac, sex_chr_dict=sex_chr_dict, faidx_dicts=faidx_dicts_cmac, chr=chr,verbose=False, log_transform=logtr)
            print(perc_X_overlap)
            plot_heatmap(counts_array=perc_X_overlap, species_list=species_list, log_data=logtr,
            filename = f"{outdir_wga}/{chr}_chr_{log_text}_alignment_coverage_heatmap_Cmac_populations.png", title = f"{chr}-Chromosome aln. coverage")

