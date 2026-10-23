"""
check for median differences in log2FC between genes sorted by different categories of gene families they belong to
"""


import parse_orthogroups as og
import parse_gff as gff
from sex_chromosomes import get_contig_names
import plot_gene_counts as data_paths
from GF_sizes import parse_orthogroups_with_gff_class

def logFC_files_dict(username="miltr339"):
    datadir = f"/Users/{username}/work/PhD_code/PhD_chapter2/data/DE_analysis"
    outdict = {
        f"{datadir}/Cmac_sex_DE_abdomen_edgeR.txt",
        f"{datadir}/Cmac_sex_DE_head_thorax_edgeR.txt",
    }
    return outdict


def plot_logFC_change(hog_sexchr_dict, parsed_gff_dict, logFC_dict, sexchr = "X", outfile = "", species = "C_maculatus"):
    """
    partition the geneIDs by their gene family assignment into groups and plot the logFC boxplots for them
    sort into: 
    * excl. A-linked
        - all genes in excl. A-linked gene families
    * X-linked (not Y):
        - all A genes in X-linked gene families
        - X genes of gene families that are X-linked in all species
        - X genes of gene families that are X-linked only in the focal species
    * x and Y linked:
        - A genes of gene families that are X and Y linked
        - X genes that are part of a gene family that also has members on the Y.
    """

    ## sort orthogroups into classes
    other_sexchr = "Y"
    if sexchr == "Y":
        other_sexchr = "X"
    hog_class_dict = {"A" : [], f"{sexchr}_all" : [], f"{sexchr}_{species}" : [], f"{sexchr}{other_sexchr}_{species}" : []}

    for hog , species_sexchr_linkage in hog_sexchr_dict.items():
        hog_species = list(species_sexchr_linkage.keys())
        # count in how many species the HOG is sex-linked
        sexchr_link_count = sum([1 for s,linkd in species_sexchr_linkage.items() if linkd[sexchr] >0 and linkd[other_sexchr]==0])
        A_link_count = sum([1 for s,linkd in species_sexchr_linkage.items() if linkd["O"] >0 and linkd[sexchr]==0 and linkd[other_sexchr]==0])
        
        if sexchr_link_count == len(hog_species):
            hog_class_dict[f"{sexchr}_all"].append(hog)
        if A_link_count == len(hog_species):
            hog_class_dict["A"].append(hog)
        if sexchr_link_count == 1 and species_sexchr_linkage[species][sexchr]>=1:
            hog_class_dict[f"{sexchr}_{species}"].append(hog)
        if species_sexchr_linkage[species][sexchr]>=1 and species_sexchr_linkage[species][other_sexchr]>=1:
            hog_class_dict[f"{sexchr}{other_sexchr}_{species}"].append(hog)
    
    for c,l in hog_class_dict.items():
        print(f"{c} : {len(l)} HOGs")



    


if __name__ == "__main__":

    username="miltr339"
    outdir = f"/Users/{username}/work/PhD_code/PhD_chapter2/data/DE_analysis"
    logfc_files = logFC_files_dict(username=username)

    of_data_dir = f"/Users/{username}/work/PhD_code/PhD_chapter2/data/orthofinder"
    orthogroups_file = f"{of_data_dir}/N0.tsv" 
    
    sex_chromosome_contigs = get_contig_names()
    annot_dict = data_paths.annotations_dict(username=username)
    
    species_order = [
        "C_maculatus",
        "C_chinensis",
        "A_obtectus",
        "B_siliquastri",
        "B_varius",
        "D_carinulata",
        "D_sublineata"
    ]

    if True:
        incl_miniprot = False

        if incl_miniprot:
            miniprot_dict = data_paths.get_miniprot_paths(username=username)
        else:
            miniprot_dict = {}

        hog_sexchr_dict, parsed_gff_dict = parse_orthogroups_with_gff_class(
            filepath=orthogroups_file, 
            annotations_dict=annot_dict,  
            sex_chr_dict=sex_chromosome_contigs, 
            miniprot_paths_dict=miniprot_dict,
            return_annotdict=True)


        for sexchr in ["X", "Y"]:
            print(f"\n\n>>>>> {sexchr} <<<<<")
            plot_logFC_change(
                hog_sexchr_dict=hog_sexchr_dict, 
                parsed_gff_dict=parsed_gff_dict, 
                logFC_dict=logfc_files,
                sexchr = sexchr, species="C_maculatus",
                outfile=f"{outdir}/logFC_boxplot_{sexchr}.png"
            )
        