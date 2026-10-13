"""
Plot gene counts in each species on the sex chromosomes
"""

import matplotlib.pyplot as plt
import parse_gff as gff
from Bio import SeqIO, Phylo, SeqUtils
from matplotlib.ticker import FuncFormatter
import numpy as np
import subprocess as sp
import pandas as pd
from sex_chromosomes import get_contig_names
import miniprot_stats_comparison as minialn
from pairwise_wga_heatmap import fasta_indices

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

def get_miniprot_paths(username="miltr339"):
    dirpath = f"/Users/{username}/work/PhD_code/PhD_chapter2/data/miniprot/"
    outdict = {
        "A_obtectus" : f"{dirpath}A_obtectus_miniprot_no_cross_no_self_hits.gff",
        "B_siliquastri" : f"{dirpath}B_siliquastri_miniprot_no_cross_no_self_hits.gff",
        "B_varius" : f"{dirpath}B_varius_miniprot_no_cross_no_self_hits.gff",
        "C_chinensis" : f"{dirpath}C_chinensis_miniprot_no_cross_no_self_hits.gff",
        "C_maculatus" : f"{dirpath}C_maculatus_miniprot_no_cross_no_self_hits.gff",
        "D_carinulata" : f"{dirpath}D_carinulata_miniprot_no_cross_no_self_hits.gff",
        "D_sublineata" : f"{dirpath}D_sublineata_miniprot_no_cross_no_self_hits.gff",
    }
    return outdict


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

def get_gene_counts_from_miniprot(miniprot_path, contig_list):
    mini_dict = minialn.miniprot_parse_alignment(miniprot_path, queryIDs=False)
    no_genes=0
    for mini_ID,aln_class in mini_dict.items():
        if aln_class.contig in contig_list:
            no_genes+=1
    return no_genes

def get_coords(tree):
    depths = tree.depths()                      # cumulative branch length from root
    if not max(depths.values()):                # no branch lengths -> cladogram
        depths = tree.depths(unit_branch_lengths=True)
    max_depth = max(depths[t] for t in tree.get_terminals())

    coords = {}
    for i, leaf in enumerate(tree.get_terminals()):
        coords[leaf] = (i, depths[leaf])

    # internal nodes: x = midpoint of children, y = own depth
    for clade in tree.get_nonterminals(order="postorder"):
        xs = [coords[c][0] for c in clade.clades]
        coords[clade] = ((min(xs) + max(xs)) / 2, depths[clade])
    return coords, max_depth



def plot_tree_manually(species_tree, lw, ax_tree = None, add_leaf_label=False, exclude_leaves=None):
    """
    plot a phylogenetic tree (newick format in a file in species_tree)
    manually so that the leaves point upwards. 
    ax_tree is the plot axis defined in fig, (ax_data, ax_tree) = plt.subplots()
    """
    tree = Phylo.read(species_tree, "newick")
    Phylo.draw_ascii(tree)

    # Prune excluded leaves before computing coordinates
    if exclude_leaves is not None:
        if isinstance(exclude_leaves, str):
            exclude_leaves = [exclude_leaves]
        for name in exclude_leaves:
            try:
                tree.prune(name)
            except ValueError:
                print(f"Warning: leaf '{name}' not found in tree, skipping")

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
                ax_tree.plot([x2, x2], [y, y2], color="black", linewidth=lw)
                # Horizontal line
                ax_tree.plot([x, x2], [y, y], color="black", linewidth=lw)

        # Draw labels on tips
        if clade.is_terminal():
            if ax_tree != None:
                if not add_leaf_label:
                    ax_tree.plot([x, x], [y, y + max_depth-y + 0.2], color="black", linewidth=lw)
                if add_leaf_label:
                    pad = 0.02 * max_depth          # small gap above the tip, scales with the tree
                    # ax_tree.text(x, y + pad, clade.name, ha="center", va="bottom", rotation=90)
                    ax_tree.plot([x, x], [y, max_depth], color="black", linestyle=":", linewidth=lw*0.75)
            leaf_names[x] = clade.name

    # Adjust and flip y-axis so tree grows upwards
    if ax_tree != None:
        ax_tree.set_ylim(ax_tree.get_ylim()[::-1])  # Invert Y axis
        ax_tree.invert_yaxis()
        ax_tree.axis("off")

    return leaf_names



def plot_gene_counts(annot_dict, species_tree, sex_chromosomes_dict, miniprot_dict = {}, chr_list=["X","Y"], filename = "only_genome_sizes_14_species.png"):
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
    lw = 4

    if type(species_tree) != list:
        # species_names = gff.make_species_order_from_tree(species_tree)
        fig, (ax_data, ax_tree) = plt.subplots(2, 1, figsize=(9, 18), gridspec_kw={'height_ratios': [1, 3]})
    
        # plot tree manually with leaves pointing upwards
        species_names_unsorted = plot_tree_manually(species_tree, lw=lw, ax_tree=ax_tree, add_leaf_label=True, exclude_leaves="C_maculatusC")
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
    if miniprot_dict != {}:
        miniprot_nos_chr = {
            chr : [get_gene_counts_from_miniprot(miniprot_path=miniprot_dict[species], contig_list=sex_chromosomes_dict[species][chr])+gene_nos_chr[chr][i] for i,species in enumerate(species_names)] for chr in chr_list
        }

    print("----------------------------------------------------")
    print(f"GENE NUMBERS\n - species\t\tchromosome: annotated + miniprot = sum\n")
    for i, species in enumerate(species_names):
        for chr in ["X","Y"]:
            print(f" - {species}\t\t{chr}: {gene_nos_chr[chr][i]} + {miniprot_nos_chr[chr][i]-gene_nos_chr[chr][i]} = {miniprot_nos_chr[chr][i]}")
    print("----------------------------------------------------")

    chr = chr_list[0]
    ax_data.plot(species_names, gene_nos_chr[chr], label = f"{chr} chromosome", color = colors[chr], linewidth = lw) 
    if miniprot_dict!={}:
        ax_data.plot(species_names, miniprot_nos_chr[chr], label = f"{chr} with miniprot", linestyle=":", color = colors[chr], linewidth = lw) 
    ylab = f"annotated genes on {chr}"
    ax_data.set_ylabel(ylab, color = colors[chr], fontsize = fs)
    ax_data.tick_params(axis ='y', labelcolor = colors[chr], labelsize = fs)  
    ymin,ymax = ax_data.get_ylim()
    ax_data.set_ylim(0,ymax)
    
    chr=chr_list[1]
    ax2 = ax_data.twinx() 
    ax2.plot(species_names, gene_nos_chr[chr], label = f"{chr} chromosome", color = colors[chr], linewidth = lw) 
    if miniprot_dict!={}:
        ax2.plot(species_names, miniprot_nos_chr[chr], label = f"{chr} with miniprot", linestyle=":", color = colors[chr], linewidth = lw) 
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


def get_window_intervals(contig_length, window_length:int) -> list[list[int]]:
    """
    get the intervals of the non-overlapping windows in a list of lists like this:
    [ [window1_start, window1_end] ,  [window2_start, window2_end] ,  ... ]
    """
    contig_end = contig_length
    contig_start = 0

    windows_list: list[list[int]] = []
    curr_window_start = contig_start
    curr_window_end = window_length

    while curr_window_end < contig_end:
        windows_list.append([curr_window_start, curr_window_end])
        curr_window_start = curr_window_end+1
        curr_window_end = curr_window_start+window_length-1
    
    windows_list.append([curr_window_start, contig_end])

    return(windows_list)




def plot_X_gene_counts_window(orthologs_path, diorhabda_annotations, sex_chromosomes, diorhabda_indices, plot_categories, winlen = 1e5, outfile_prefix=None):
    """
    Plot gene numbers in winlen bp long windows in two categories: X in Diorhabda, A everywhere else, and X everywhere
    """
    plt.rcParams['text.usetex'] = True # use \\textit{{{}}} for species names
    plt.rcParams['text.latex.preamble'] = r'\usepackage{sfmath} \renewcommand{\familydefault}{\sfdefault}'
    plt.rcParams['font.family'] = 'sans-serif'

    colors = {
        "X in bruchini" : "#776274", # dusty lavender
        "A in bruchini" : "#006D77", # stormy teal
    }
    fs = 30 # set font size
    lw = 3

    orthologs_df = pd.read_csv(orthologs_path)

    linkage_df_dict = {}
    for chr_linkage, linkage_dict in plot_categories.items():
        mod_df = orthologs_df
        for h,chr_ in linkage_dict.items():
            mod_df = mod_df[mod_df[h]==chr_]
        linkage_df_dict[chr_linkage] = mod_df

        chr_linkage_ = chr_linkage.replace(" ","_")
        orthologs_path_ = orthologs_path.replace(".csv", f"_{chr_linkage_}.csv")
        ## save orthogroups to file for dNdS without including the columns that have the A or X category
        mod_df.to_csv(path_or_buf=orthologs_path_, index=False , columns=[c for c in mod_df.columns if ":chr" not in c])

        print(f"  - {mod_df.shape[0]} orthologs in '{chr_linkage}'")
    
    for diorhabda_species, d_annotation in diorhabda_annotations.items():
        
        filename = outfile_prefix.replace(".png", f"_{diorhabda_species}.png")
        fig, ax = plt.subplots(1, 1, figsize=(13, 8))
        X_chr = sex_chromosomes[diorhabda_species]["X"][0]

        print(f"parse annotation: {d_annotation}\nonly including the X chromosome '{X_chr}'")
        annot_df = gff.parse_gff3_general(d_annotation, keep_feature_category=gff.FeatureCategory.Transcript, contig_list=[X_chr], verbose=False)
        assert len(annot_df)>0

        ## get window intervals, compute window edges
        window_intervals = get_window_intervals(contig_length=diorhabda_indices[diorhabda_species][X_chr], window_length=winlen)
        
        for plot_cat, data_df in linkage_df_dict.items():
            legend_label = plot_cat.replace("bruchini", "\\textit{{{Bruchini}}}")
            transcriptIDs = data_df[diorhabda_species].values.tolist()
            x_coords = []
            y_coords = []
            for window_interval in window_intervals:
                window_start,window_stop = window_interval
                x_coords.append(window_start)
                x_coords.append(window_stop)
                genecount_window = 0
                ## I count a gene as present in a window by the START of the feature annotation! if a gene goes across a window border it is counted in the earlier window
                for geneID in transcriptIDs:
                    if annot_df[geneID].start < window_stop and annot_df[geneID].start > window_start:
                        genecount_window +=1
                # append twice to match x coords
                y_coords.append(genecount_window)
                y_coords.append(genecount_window)
            
            ax.plot(x_coords, y_coords, label = legend_label, color = colors[plot_cat], linewidth = lw) 
            plt.fill_between(x=x_coords, y1=y_coords, color=colors[plot_cat], alpha = 0.4)

        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: '' if x <= 0 else f'{x / 1e6:.0f} Mb'))
        ax.set_ylabel(f"gene count per window", fontsize = fs)
        ax.tick_params(axis='y', labelsize=fs)
        ax.tick_params(axis='x', labelsize=fs)
        ax.set_xlim(0,diorhabda_indices[diorhabda_species][X_chr])
        ax.set_ylim(0,29)

        species = diorhabda_species.replace("_", ". ")
        plt.title(f"Gene counts on \\textit{{{species}}} X chromosome in {winlen / 1e3:.0f}kb windows", fontsize = fs)
        plt.legend(fontsize=fs)
        plt.tight_layout()
        plt.savefig(filename, dpi = 300, transparent = True)# , bbox_inches='tight')

        print(f"{filename}")



if __name__=="__main__":

    username="miltr339"
    sex_chromosome_contigs = get_contig_names()
    annot_dict = annotations_dict(username=username)
    faidx_dict,_ = fasta_indices(username=username)
    miniprot_dict = get_miniprot_paths(username=username)
    species_order = [
"C_maculatus",
"C_chinensis",
"A_obtectus",
"B_siliquastri",
"B_varius",
"D_carinulata",
"D_sublineata"]
    species_tree_path=f"/Users/{username}/work/PhD_code/PhD_chapter2/data/orthofinder/SpeciesTree_rooted_node_labels.txt"

    if False:
        data_dir = f"/Users/{username}/work/PhD_code/PhD_chapter2/data/"
        plot_gene_counts(
            annot_dict=annot_dict, 
            # species_tree=species_order, 
            species_tree=species_tree_path, 
            sex_chromosomes_dict=sex_chromosome_contigs, 
            chr_list=["X","Y"], 
            miniprot_dict=miniprot_dict,
            filename = f"{data_dir}sex_chromosome_gene_counts_tree.png")

    if True:
        # Plot gene counts in windows for Diorhabda X chromosomes
        diorhabda = ["D_carinulata","D_sublineata"]
        bruchini = ["C_maculatus","C_chinensis","A_obtectus","B_siliquastri","B_varius"]
        orthologs_path=f"/Users/{username}/work/PhD_code/PhD_chapter2/data/BRH_orthologs/chrysomelidae_orthologs.csv"
        plot_cats = {
            "X in bruchini" : {f"{s}:chr" : "X" for s in species_order},
            "A in bruchini" : {f"{s}:chr" : "X" if s in diorhabda else "A" for s in species_order}
        }
        print(plot_cats)
        plot_X_gene_counts_window(
            orthologs_path=orthologs_path, 
            diorhabda_annotations ={d : annot_dict[d] for d in diorhabda}, 
            diorhabda_indices = {d : faidx_dict[d] for d in diorhabda},
            sex_chromosomes = sex_chromosome_contigs, 
            plot_categories = plot_cats, winlen = 75e4, 
            outfile_prefix=f"/Users/{username}/work/PhD_code/PhD_chapter2/data/BRH_orthologs/XA_gene_density.png")