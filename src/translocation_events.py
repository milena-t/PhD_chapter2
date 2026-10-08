import linkage_groups as lg
import GF_sizes as OGs
import linkage_groups as lg
import parse_orthogroups as og
import plot_gene_counts as data_paths
from sex_chromosomes import get_contig_names 

"""
Categorize gene family translocation and duplication events.
"""



if __name__ == "__main__":
    
    username="miltr339"

    data_dir = f"/Users/{username}/work/PhD_code/PhD_chapter2/data/orthofinder"
    orthogroups_file = f"{data_dir}/N0.tsv" 
    unassigned_genes_path = f"{data_dir}/unassigned_genes.tsv" 
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

    orthogorups = og.parse_orthogroups_class(
        filepath=orthogroups_file,
        annot_species=species_order,
        annotations_dict=annot_dict,
        unassigned_genes_path=unassigned_genes_path,
        miniprot_paths_dict=miniprot_dict)
    
    print(orthogorups["N0.HOG0000064"])