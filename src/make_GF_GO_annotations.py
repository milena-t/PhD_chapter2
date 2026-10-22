"""
I will make gene-family level GO annotations for each species based on the transcript-level annotations from eggnog
I will summarize all GO terms for each gene family and then just leave a list of unique ones
"""

import parse_orthogroups as og
from sex_chromosomes import get_contig_names
import plot_gene_counts as data_paths
import pandas as pd
from collections import defaultdict

def eggnog_paths(username="miltr339"):
    datadir = f"/Users/{username}/work/chapter2/eggnog/"
    outdict = {
        "A_obtectus" : f"{datadir}A_obtectus_eggnog_diamond.emapper.annotations",
        "B_siliquastri" : f"{datadir}B_siliquastri_eggnog_diamond.emapper.annotations",
        "B_varius" : f"{datadir}B_varius_eggnog_diamond.emapper.annotations",
        "C_chinensis" : f"{datadir}C_chinensis_eggnog_diamond.emapper.annotations",
        "C_maculatus" : f"{datadir}C_maculatus_eggnog_diamond.emapper.annotations",
        "D_carinulata" : f"{datadir}D_carinulata_eggnog_diamond.emapper.annotations",
        "D_sublineata" : f"{datadir}D_sublineata_eggnog_diamond.emapper.annotations",
        "T_castaneum" : f"{datadir}T_castaneum_eggnog_diamond.emapper.annotations",
    }
    return outdict

def make_gf_eggnog(eggnog_path, trID_hog_lookup, outfile_name):
    """
    read the eggnog output into a df, assign 
    """
    header_list = ["#query","seed_ortholog","evalue","score","eggNOG_OGs","max_annot_lvl","COG_category","Description","Preferred_name","GOs","EC","KEGG_ko","KEGG_Pathway","KEGG_Module","KEGG_Reaction","KEGG_rclass","BRITE","KEGG_TC","CAZy","BiGG_Reaction","PFAMs"]
    GO_header = "GOs"
    df = pd.read_csv(eggnog_path, sep="\t", comment="#", names=header_list, index_col=0)

    hog_eggnog_dict = defaultdict(set)
    for trid, hog in trID_hog_lookup.items():
        try:
            row = df.loc[trid]
        except:
            continue
        GOs = row[GO_header]
        if GOs == "-":
            continue
        GOs = set(GOs.split(","))
        hog_eggnog_dict.setdefault(hog, set()).update(GOs)

    print(f"{len(trID_hog_lookup)} transcript IDs, {len(set(trID_hog_lookup.values()))} orthogroups, {len(hog_eggnog_dict)} of which have functional annotations")
    hog_eggnog_df = pd.DataFrame({hog : {"GOs" : ",".join(gos)} for hog,gos in hog_eggnog_dict.items()}).T
    hog_eggnog_df.index.name = "#query"

    # add empty cols for other functional categories
    header_list.remove("#query")
    header_list_all = [h for h in header_list]
    header_list.remove(GO_header)
    for header in header_list:
        hog_eggnog_df[header] = "-"
    # reorder columns to match eggnog output
    hog_eggnog_df = hog_eggnog_df[header_list_all]

    # save as list
    hog_eggnog_df.to_csv(outfile_name, sep="\t")




def get_orthogroups_lookup_dict(orthogroups_path, species):
    """
    make a lookup dict for each species to check OG membership of each transcript
    """
    HOGs_dict = og.parse_orthogroups_dict(filepath=orthogroups_path, species_list=[species])
    HOGs_dict = {hog : member[species] for hog,member in HOGs_dict.items()}

    outdict = {}
    for hog, trIDs_str in HOGs_dict.items():
        if type(trIDs_str)==float:
            continue ## true if the variable is nan
        outdict = outdict | {trid : hog for trid in trIDs_str.split(", ")}
            
    return outdict


if __name__ == "__main__":
    
    username="miltr339"

    data_dir = f"/Users/{username}/work/PhD_code/PhD_chapter2/data/orthofinder"
    orthogroups_file = f"{data_dir}/N0.tsv" 
    outdir = f"/Users/{username}/work/PhD_code/PhD_chapter2/data/GO_enrichment/"
    unassigned_genes_path = f"{data_dir}/unassigned_genes.tsv" 
    sex_chromosome_contigs = get_contig_names()
    annot_dict = data_paths.annotations_dict(username=username)
    miniprot_dict = data_paths.get_miniprot_paths(username=username)
    eggnog_dict = eggnog_paths(username=username)
    species_order = [
        "C_maculatus",
        "C_chinensis",
        "A_obtectus",
        "B_siliquastri",
        "B_varius",
        "D_carinulata",
        "D_sublineata"
    ]


    for species in species_order:
        print(f"\n---<>--> {species}")
        trIDs_lookup = get_orthogroups_lookup_dict(orthogroups_path = orthogroups_file, species=species)
        make_gf_eggnog(eggnog_path=eggnog_dict[species], trID_hog_lookup=trIDs_lookup, outfile_name=f"{outdir}{species}_HOG_GO_annotations.tsv")
        