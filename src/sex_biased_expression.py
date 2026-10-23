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


def plot_logFC_change(hog_sexchr_dict, parsed_gff_dict, logFC_dict, N0_file, species_sexchr_dict, sexchr = "X", outfile = "", species = "C_maculatus",):
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
        if species not in hog_species:
            continue
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
    print(f"---------------------------------")

    ## sort genes from orthogroup classes
    annot_placed_dict = og.parse_orthogroups_dict(filepath=N0_file, species_list=[species], OG_header="HOG")
    gene_class_dict = {
        "A_in_A_gf" : [],
        
        f"A_in_{sexchr}_all_gf" : [],
        f"X_in_{sexchr}_all_gf" : [],        ## always genes on X !! sex-bias in Y is not very informative
        
        f"A_in_{sexchr}_{species}_gf" : [],
        f"X_in_{sexchr}_{species}_gf" : [],  ## always genes on X !! sex-bias in Y is not very informative

        f"A_in_{sexchr}{other_sexchr}_all_gf" : [],
        f"X_in_{sexchr}{other_sexchr}_{species}_gf" : [],  ## always genes on X !! sex-bias in Y is not very informative
    }

    for hog_category, hog_list in hog_class_dict.items():
        for hog in hog_list:
            for sp,tr_string in annot_placed_dict[hog].items():
                assert sp == species
                tr_list = tr_string.split(", ")
                for trID in tr_list:
                    tr_contig = "A"
                    for sc, clist in species_sexchr_dict.items():
                        if parsed_gff_dict[trID].contig in clist:
                            tr_contig=sc

                    if tr_contig=="Y":
                        continue
                    gID = parsed_gff_dict[trID].parent_id

                    if hog_category == "A" and tr_contig == "A":
                        gene_class_dict["A_in_A_gf"].append(gID)

                    if hog_category == f"{sexchr}_all" and tr_contig == "A":
                        gene_class_dict[f"A_in_{sexchr}_all_gf"].append(gID)
                    if hog_category == f"{sexchr}_all" and tr_contig == "X":
                        gene_class_dict[f"X_in_{sexchr}_all_gf"].append(gID) 
                    if hog_category == f"{sexchr}_{species}" and tr_contig == "A":
                        gene_class_dict[f"A_in_{sexchr}_{species}_gf"].append(gID)
                    if hog_category == f"{sexchr}_{species}" and tr_contig == "X":
                        gene_class_dict[f"X_in_{sexchr}_{species}_gf"].append(gID) 
                    if hog_category ==  f"{sexchr}{other_sexchr}_{species}" and tr_contig == "A":
                        gene_class_dict[f"A_in_{sexchr}{other_sexchr}_all_gf"].append(gID)
                    if hog_category ==  f"{sexchr}{other_sexchr}_{species}" and tr_contig == "X":
                        gene_class_dict[f"X_in_{sexchr}{other_sexchr}_{species}_gf"].append(gID) 

    for c,l in gene_class_dict.items():
        print(f"{c} : {len(l)} genes, {l[:10]}...")
    print(f"---------------------------------")

    



    


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

        hog_sexchr_dict,_ = parse_orthogroups_with_gff_class(
            filepath=orthogroups_file, 
            annotations_dict=annot_dict,  
            sex_chr_dict=sex_chromosome_contigs, 
            miniprot_paths_dict=miniprot_dict,
            return_annotdict=True)

        focal_species ="C_maculatus"
        species_full_GFFdict = gff.parse_gff3_general(filepath=annot_dict[focal_species], verbose=False)
        for sexchr in ["X", "Y"]:
            print(f"\n\n>>>>>>>>>>>>>>> {sexchr} <<<<<<<<<<<<<<<")
            plot_logFC_change(
                hog_sexchr_dict=hog_sexchr_dict, 
                parsed_gff_dict=species_full_GFFdict, 
                logFC_dict=logfc_files,
                species_sexchr_dict = sex_chromosome_contigs[focal_species],
                sexchr = sexchr, species=focal_species,
                N0_file = orthogroups_file,
                outfile=f"{outdir}/logFC_boxplot_{sexchr}.png"
            )
        