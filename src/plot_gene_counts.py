"""
Plot gene counts in each species on the sex chromosomes
"""

import matplotlib.pyplot as plt
import parse_gff as gff
from Bio import SeqIO, Phylo, SeqUtils
import numpy as np
import subprocess as sp
from sex_chromosomes import get_contig_names

def annotations_dict(username="miltr339"):
    dirname = f"/Users/{username}/work/chapter2/native_annotations"
    # /proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/native_annotations
    out_dict = {
        "A_obtectus" : f"{dirname}/A_obtectus.gff",
        "B_siliquastri" : f"{dirname}/B_siliquastri.gff",
        "B_varius" : f"{dirname}/B_varius.gff",
        "C_chinensis" : f"{dirname}/C_chinensis.gff",
        "C_maculatus" : f"{dirname}/C_maculatus.gff",
        "D_carinulata" : f"{dirname}/D_carinulata.gff",
        "D_sublineata" : f"{dirname}/D_sublineata.gff",
    }   
    return out_dict


def get_gene_conuts_from_annot(annot_path, contig_list):
    
    try:
        annot_dict = gff.parse_gff3_by_contig(annot_path)
    except:
        annot_dict = gff.parse_gff3_by_contig(annot_path, gtf=True)

    protein_num_list = []
    no_genes = 0
    for contig in contig_list:
        if contig in annot_dict:
            protein_num_list.append(len(annot_dict[contig]))
        else:
            no_genes+=1
    num_proteins = sum(protein_num_list)
    # print(f"{no_genes} contigs have no annotated genes")
    return num_proteins


def get_coords(tree):
    """
    Traverse tree and assign x/y positions for each clade.
    The tree orientation is intended with leaves pointing upwards in a normal x/y coordinate system
    """
    coords = {}
    x = 0
    max_depth = 0 # get position of the deepest leaf so that you can align all the leaf labels later
    leaf_names = []
    
    # recursion through every node of the tree
    def assign(clade, depth):
        nonlocal x
        nonlocal max_depth
        
        if depth>max_depth:
            max_depth = depth
        if clade.is_terminal():
            coords[clade] = (x, depth) # depth is the distance from the root
            x += 1 # orientation of the node along the x-axis (left/right orientation)
        else:
            for child in clade.clades:
                assign(child, depth + clade.branch_length if clade.branch_length else depth)
            # if the node has multiple children, determine the x coordinates of all the children and take the mean to center the line in the middle above the children
            child_coords = [coords[c][0] for c in clade.clades]
            coords[clade] = (sum(child_coords) / len(child_coords), depth)

    assign(tree.root, 0)

    return coords, max_depth



def plot_tree_manually(species_tree, ax_tree = None, add_leaf_label=False):
    """
    plot a phylogenetic tree (newick format in a file in species_tree)
    manually so that the leaves point upwards. 
    ax_tree is the plot axis defined in fig, (ax_data, ax_tree) = plt.subplots()
    """
    tree = Phylo.read(species_tree, "newick")
    tree.ladderize()
    coords, max_depth = get_coords(tree)
    leaf_names = {}
    # Draw manually
    for clade in tree.find_clades(order="level"):
        x, y = coords[clade]

        # Draw a horizontal line from this node to its children
        for child in clade.clades:
            x2, y2 = coords[child]

            if ax_tree != None:
                # Vertical line
                ax_tree.plot([x2, x2], [y, y2], color="black")
                # Horizontal line
                ax_tree.plot([x, x2], [y, y], color="black")

        # Draw labels on tips
        if clade.is_terminal():
            if ax_tree != None:
                ax_tree.plot([x, x], [y, y + max_depth-y + 0.2], color="black")
                if add_leaf_label:
                    ax_tree.text(x, y + max_depth-y + 0.25, clade.name, ha="center", va="bottom", rotation=90)
            leaf_names[x] = clade.name

    # Adjust and flip y-axis so tree grows upwards
    if ax_tree != None:
        ax_tree.set_ylim(ax_tree.get_ylim()[::-1])  # Invert Y axis
        ax_tree.invert_yaxis()
        ax_tree.axis("off")

    return leaf_names



def plot_gene_counts(annot_dict, species_tree, sex_chromosomes_dict, chr_list=["X","Y"], filename = "only_genome_sizes_14_species.png"):
    """
    plot gene counts from annotations (or proteinfasta, but preferably annotation), with a species tree on the x-axis
    """
    plt.rcParams['text.usetex'] = True # use \\textit{{{}}} for species names
    plt.rcParams['text.latex.preamble'] = r'\usepackage{sfmath} \renewcommand{\familydefault}{\sfdefault}'
    plt.rcParams['font.family'] = 'sans-serif'

    colors = {
        "X" : "#B12E17", # red
        "Y" : "#768DC4", # blue
    }
    
    fs = 30 # set font size

    if type(species_tree) != list:
        species_names = gff.make_species_order_from_tree(species_tree)
        fig, (ax_data, ax_tree) = plt.subplots(2, 1, figsize=(10, 15), gridspec_kw={'height_ratios': [1, 2]})
    
        # plot tree manually with leaves pointing upwards
        species_names_unsorted = plot_tree_manually(species_tree, ax_tree)
        # get species order from plotted tree
        species_coords_sorted = sorted(list(species_names_unsorted.keys()))
        species_names = [species_names_unsorted[species_coord] for species_coord in species_coords_sorted]
    else:
        fig, ax_data = plt.subplots(1, 1, figsize=(13, 8))
        species_names = species_tree

    # get a list of lists with [native, orthoDB] number of gene families per species

    legend_labels = []
    ncol_legend = 0

    gene_nos_chr = {
        chr : [get_gene_conuts_from_annot(annot_dict[species], contig_list=sex_chromosomes_dict[species][chr]) for species in species_names] for chr in chr_list
    }

    chr = chr_list[0]
    ax_data.plot(species_names, gene_nos_chr[chr], label = f"{chr} chromosome", color = colors[chr], linewidth = 4) # red
    ylab = f"annotated genes on {chr}"
    ax_data.set_ylabel(ylab, color = colors[chr], fontsize = fs)
    ax_data.tick_params(axis ='y', labelcolor = colors[chr], labelsize = fs)  
    ymin,ymax = ax_data.get_ylim()
    ax_data.set_ylim(0,ymax)
    
    chr=chr_list[1]
    ax2 = ax_data.twinx() 
    ax2.plot(species_names, gene_nos_chr[chr], label = f"{chr} chromosome", color = colors[chr], linewidth = 4) # red
    ylab = f"annotated genes on {chr}"
    ax2.set_ylabel(ylab, color = colors[chr], fontsize = fs)
    ax2.tick_params(axis ='y', labelcolor = colors[chr], labelsize = fs)  
    ymin,ymax = ax2.get_ylim()
    ax2.set_ylim(0,ymax*1.25)

    species_axis_labels = [species.replace("_", ". ") for species in species_names]
    ax_data.set_xticklabels([f"\\textit{{{species}}}" for species in species_axis_labels], rotation=90, fontsize=fs)
    
    # yhandles, ylabels = ax2.get_legend_handles_labels()
    # xhandles, xlabels = ax_data.get_legend_handles_labels()
    # legend = ax_data.legend(handles = [yhandles,xhandles], labels=[ylabels,xlabels], fontsize = fs*0.8)

    # set grid only for X axis ticks 
    ax_data.grid(True)
    ax_data.yaxis.grid(False)
    
    ax_data.tick_params(axis='y', labelsize=fs)

    plt.tight_layout()

    plt.savefig(filename, dpi = 300, transparent = True)# , bbox_inches='tight')
    print("Figure saved in the current working directory directory as: "+filename)


if __name__=="__main__":

    username="miltr339"
    sex_chromosome_contigs = get_contig_names()
    annot_dict = annotations_dict(username=username)
    species_order = [
"C_maculatus",
"C_chinensis",
"A_obtectus",
"B_siliquastri",
"B_varius",
"D_carinulata",
"D_sublineata"]
    data_dir = f"/Users/{username}/work/PhD_code/PhD_chapter2/data/"
    plot_gene_counts(annot_dict, species_tree=species_order, sex_chromosomes_dict=sex_chromosome_contigs, chr_list=["X","Y"], filename = f"{data_dir}sex_chromosome_gene_counts.png")