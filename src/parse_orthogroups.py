# check the position of the significant orthogroups in each annotation:
#  - do they look like tandem duplications or are they dispersed?
#  - Are the n basepairs around the gene enriched for a certain kind of TE compared to the rest of the genome?
#    (pick n like 1e5? pretty short probably)
#  - Are they intron-less?

import upsetplot
import numpy as np
import pandas as pd
from statistics import mean
import matplotlib.pyplot as plt

import parse_gff as gff
from sex_chromosomes import get_contig_names
import miniprot_stats_comparison as minialn
import plot_gene_counts as data_paths



def get_orthogroup_sizes(orthogroup_dict, q = 0):
    """
    returns a dict of all orthogroups sizes. Do NOT confuse it with gene family sizes! (get_GF_sizes function instead)
    if a percentile q is specified other than 0, it returns only orthogroups whose size is > than the pth percentile of the size distribution
    """
    
    OG_sizes_dict = {}
    sizes = []
    for OG_id, species_dict in orthogroup_dict.items():
        size = 0
        for transcripts_list in species_dict.values():
            size += len(transcripts_list)
        
        OG_sizes_dict[OG_id] = size
        sizes.append(size)
    
    if q == 0:
        return OG_sizes_dict
    else:
        OG_sizes_filtered = {}
        sizes = np.array(sizes)
        percentile_size = np.percentile(sizes, q = q)
        for OG_id, size in OG_sizes_dict.items():
            if size > percentile_size:
                OG_sizes_filtered[OG_id] = size

        return OG_sizes_filtered


def parse_orthogroups_dict(filepath, species_list):
    orthogroups_df = pd.read_csv(filepath, sep="\t")
    headers = set(orthogroups_df.columns)
    annot_species = set(species_list)
    headers_keep = list(annot_species & headers)
    headers = ["HOG"] + headers_keep
    
    df = orthogroups_df[headers]
    df = df.set_index("HOG", drop=True)
    df_dict = df.to_dict(orient="index")
    
    return df_dict

def get_GF_sizes(orthogroups_dict):
    """
    modify the dictionary from the parse_orthogroups function
    
    {
        orthogroup : {
            species1 : num,
            species2 : num,
            ...
        }
        orthogroup2 : {
            species1 : num,
            species2 : num,
            ...
        }
    }
    """
    orthogroups_modified = {}
    for orthogroup_id , species_dict in orthogroups_dict.items():
        species_dict_nums = { species : len(transcripts_list) for species, transcripts_list in species_dict.items() if transcripts_list!=['']}
        orthogroups_modified[orthogroup_id] = species_dict_nums

    return orthogroups_modified

            
def get_mean_GF_size(orthogroups_dict:dict, species:str):
    GF_sizes = []
    for gene_families in orthogroups_dict.values():
        try:
            num_families = len(gene_families[species])
        except:
            num_families = 0
        GF_sizes.append(num_families)
    
    return mean(GF_sizes)
    
    
def get_species_in_OG_dict(OG_dict:dict) -> list:
    """
    get the list of all species included in the orthofinder dict
    """
    all_species = []
    for OG_id, GF_dict in OG_dict.items():
        species = list(GF_dict.keys())
        all_species.extend(species)
    all_species = list(set(all_species))
    return(all_species)



def parse_orthogroups_with_gff_class(filepath, annotations_dict, sex_chr_dict, miniprot_paths_dict = {}):
    """
    Get a dictionary of all the orthogroups with counts for the sex chromosomes
    dict_out = { 
        HOG_ID : {
            species_1 : { "O" : n, "X" : n, "Y" : n },
            species_2 : { "O" : m, "X" : m, "Y" : m },
            ...
        },
        ...
    }
    """

    annot_species = list(annotations_dict.keys())
    annot_HOG_dict = parse_orthogroups_dict(filepath=filepath, species_list=annot_species)

    annot_gff_dict = { species : gff.parse_gff3_general(annotations_dict[species], keep_feature_category=gff.FeatureCategory.Transcript, verbose=False) for species in annot_species}
    
    hog_sexchr_dict = {}

    if miniprot_paths_dict != {}:
        miniprot_dict = { species : minialn.miniprot_parse_alignment(miniprot_paths_dict[species]) for species in annot_species}

        for HOG_id , OG_dict in annot_HOG_dict.items():
            sexchr_counts_dict = {}
            # print(f" --- {HOG_id} --- ")
            for species, geneIDs in OG_dict.items():
                sex_chr_counts = {"O" : 0 , "X" : 0, "Y" : 0} # O for other instead of A for autosomes, because 'O' is any non-X/Y contig, also unplaced scaffolds
                if isinstance(geneIDs, str): ## if not NaN
                    geneIDs = geneIDs.strip().split(", ")
                    # print(f"   - {species} : {len(geneIDs)}")
                    for geneID in geneIDs:
                        gene_contig = annot_gff_dict[species][geneID].contig 
                        if gene_contig in sex_chr_dict[species]["X"]:
                            sex_chr_counts["X"] += 1
                        elif gene_contig in sex_chr_dict[species]["Y"]:
                            sex_chr_counts["Y"] += 1
                        else:
                            sex_chr_counts["O"] += 1
                        if geneID in miniprot_dict:
                            
                else:
                    # print(f"   - {species} : 0")
                    pass
                sexchr_counts_dict[species] = sex_chr_counts
            hog_sexchr_dict[HOG_id] = sexchr_counts_dict
    else:
        for HOG_id , OG_dict in annot_HOG_dict.items():
            sexchr_counts_dict = {}
            # print(f" --- {HOG_id} --- ")
            for species, geneIDs in OG_dict.items():
                sex_chr_counts = {"O" : 0 , "X" : 0, "Y" : 0} # O for other instead of A for autosomes, because 'O' is any non-X/Y contig, also unplaced scaffolds
                if isinstance(geneIDs, str): ## if not NaN
                    geneIDs = geneIDs.strip().split(", ")
                    # print(f"   - {species} : {len(geneIDs)}")
                    for geneID in geneIDs:
                        gene_contig = annot_gff_dict[species][geneID].contig 
                        if gene_contig in sex_chr_dict[species]["X"]:
                            sex_chr_counts["X"] += 1
                        elif gene_contig in sex_chr_dict[species]["Y"]:
                            sex_chr_counts["Y"] += 1
                        else:
                            sex_chr_counts["O"] += 1
                else:
                    # print(f"   - {species} : 0")
                    pass
                sexchr_counts_dict[species] = sex_chr_counts
            hog_sexchr_dict[HOG_id] = sexchr_counts_dict
    
    return(hog_sexchr_dict)



def upset_HOG_sex_chromosomes(hog_sexchr_dict, chr = "Y", plot_filename = "upsetplot.png", min_intersection_size = 0):

    hog_sexchr_unnested = []
    autosome_excl = 0
    chr_linked = 0
    for HOG_id, sexchr_species_dict in hog_sexchr_dict.items():
        species_incl = []
        for species, sexchr_dict in sexchr_species_dict.items():
            if sexchr_dict[chr] >0:
                species_ = species.replace("_", ". ")
                species_incl.append(f"\\textit{{{species_}}}")
        
        # only plot ones with at least one {chr}-linked member
        if species_incl !=[]:
            hog_sexchr_unnested.append(species_incl)
            chr_linked += 1
        else:
            autosome_excl+=1
        
    print(f"{autosome_excl} orthogroups autosome-exclusive, {chr_linked} are {chr}-linked in at least one species")


    plt.rcParams['text.usetex'] = True # use \textit{} for species names
    plt.rcParams['text.latex.preamble'] = r'\usepackage{sfmath} \renewcommand{\familydefault}{\sfdefault}'
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.size'] = 16
    fig = plt.figure()

    data = upsetplot.from_memberships(hog_sexchr_unnested)
    upsetplot.UpSet(data, subset_size="count", sort_by="cardinality", sort_categories_by="input", show_counts=True, min_subset_size=min_intersection_size).plot(fig=fig)

    if min_intersection_size>0:
        plot_title = f"Orthogroup presence on {chr} (intersection size $>$ {min_intersection_size})"
        fig.suptitle(plot_title, y=0.97)
    else:
        plot_title = f"Orthogroup presence on {chr}"
        fig.suptitle(plot_title)
    plt.tight_layout()
    plt.savefig(plot_filename, dpi = 300, transparent = True, bbox_inches="tight")
    print(f"plot saved in current working directory as: {plot_filename}")
    plt.clf()
    plt.cla()
    plt.close()


if __name__ == "__main__":
    
    username="miltr339"

    data_dir = f"/Users/{username}/work/PhD_code/PhD_chapter2/data/orthofinder"
    orthogroups_file = f"{data_dir}/N0.tsv" 
    sex_chromosome_contigs = get_contig_names()
    annot_dict = data_paths.annotations_dict(username=username)
    miniprot_dict = data_paths.get_miniprot_paths(username=username)
    species_order = [
        "C_maculatus",
        "C_chinensis",
        "A_obtectus",
        "B_siliquastri",
        "B_varius",
        "D_carinulata",
        "D_sublineata"
    ]
    hog_sexchr_dict = parse_orthogroups_with_gff_class(filepath=orthogroups_file, annotations_dict=annot_dict, sex_chr_dict=sex_chromosome_contigs)

    min_intersection_size = {"X" : 20, "Y" : 0}
    for chr in ["Y","X"]:
        upset_HOG_sex_chromosomes(hog_sexchr_dict=hog_sexchr_dict, chr=chr, plot_filename=f"{data_dir}/orthogroup_presence_{chr}_upsetplot.png", min_intersection_size=min_intersection_size[chr])