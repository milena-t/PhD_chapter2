import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from statistics import mean
from scipy.stats import sem
from collections import defaultdict

import linkage_groups as lg
import GF_sizes as OGs
import parse_gff as gff
import linkage_groups as lg
import parse_orthogroups as og
import plot_gene_counts as data_paths
from GF_sizes import genome_sizes

"""
Categorize gene family translocation and duplication events.
"""

def GS_vs_LG_number_correlations(orthogroups, species_order, filename_prefix, max_GF_size = 0):
    """
    plot correlation of GF size and number of LG the GF is on, split by X and Y linkage vs exclusively autosome linked
    """

    plt.rcParams['text.usetex'] = True # use \\textit{{{}}} for species names
    plt.rcParams['text.latex.preamble'] = r'\usepackage{sfmath} \renewcommand{\familydefault}{\sfdefault}'
    plt.rcParams['font.family'] = 'sans-serif'

    fs = 15 # font size
    ps = 10
    lw = 2
    aspect_ratio = 18 / 14 # height / width
    height_pixels = 1400  # Height in pixels
    dpi = 300
    width_pixels = int(height_pixels * aspect_ratio)  # Width in pixels

    colors_dict = {
        "A" : "#748B7A",
        "X" : "#BD351E", # red
        "Y" : "#5E79BD", # blue
        "XY" : "#84569E",
        "Y_edge" : "#374C6E", # dusk blue darker
        "X_edge" : "#771C2C", # dark amaranth
        "A_edge" : "#376F59", # deep teal
        "XY_edge" : "#59376D", # velvet purple
        "vline" : '#9499A5', # cool steel
        "vtext" : "#5B6378", #blue slate
    }
    legend_label = {
        "A" : "Autosomes",
        "XY" : "X and Y linked",
        "X" : "X linked",
        "Y" : "Y linked"
    }

    for species in species_order:
        # print(f"\n===================== {species} =====================")
        A_nums = {}
        X_nums = {} # GF sizes as integer keys, val is list of number of LGs
        Y_nums = {}
        XY_nums = {}

        species_=species.replace("_", ". ")
        if max_GF_size ==0:
            fig_title = f"\\textit{{{species_}}}"
        else:
            fig_title = f"\\textit{{{species_}}} (max GF size {max_GF_size})"

        for ogID, orthogroup in orthogroups.items():
            GF_size = orthogroup.GF_size(species)
            linkage_groups = orthogroup.linkagegroups(species)
            LG_number = len(linkage_groups)
            if lg.LinkageGroup.LGX in linkage_groups and lg.LinkageGroup.LGY not in linkage_groups:
                X_nums.setdefault(GF_size, []).append(LG_number)
            elif lg.LinkageGroup.LGX not in linkage_groups and lg.LinkageGroup.LGY in linkage_groups:
                Y_nums.setdefault(GF_size, []).append(LG_number)
            elif lg.LinkageGroup.LGX in linkage_groups and lg.LinkageGroup.LGY in linkage_groups:
                XY_nums.setdefault(GF_size, []).append(LG_number)
            else:
                A_nums.setdefault(GF_size, []).append(LG_number)
        
        data_dict = { # sort by key (GF size) ascending
            "A" : dict(sorted(A_nums.items())),
            "X" : dict(sorted(X_nums.items())),
            "Y" : dict(sorted(Y_nums.items())),
            "XY" : dict(sorted(XY_nums.items())),
        }
        fig, ax = plt.subplots(1,1,figsize=(width_pixels/dpi, height_pixels/dpi))

        for c,data_nums in data_dict.items():
            gf_sizes = list(data_nums.keys())
            lg_means = [0.0 for i in gf_sizes]
            lg_errs = [0.0 for i in gf_sizes]

            for i,gf_size in enumerate(gf_sizes):
                lg_means[i] = mean(data_nums[gf_size])
                lg_errs[i] = sem(data_nums[gf_size])

            ax.errorbar(gf_sizes, lg_means, yerr = lg_errs, color=colors_dict[c], linewidth =lw*0.5, marker = ".", markersize=ps, linestyle = ":", label = legend_label[c])
        
        ylab = f"Number of linkage groups"
        ax.set_ylabel(ylab, fontsize = fs)
        ax.set_xlabel("Gene family size", fontsize = fs)
        ax.tick_params(axis ='y', labelsize = fs)
        ax.tick_params(axis ='x', labelsize = fs)
        
        if max_GF_size>0:
            xmin,xmax = ax.get_xlim()
            ax.set_xlim(-1, max_GF_size)
        
        plt.title(f"{fig_title}", fontsize = fs*1.25)
        ax.legend(fontsize=fs*0.75, loc="upper left")
        plt.tight_layout()

        outfile_annot = filename_prefix.replace(".png", f"_{species}.png")
        plt.savefig(outfile_annot, dpi = 300, transparent = True)# , bbox_inches='tight')
        print("Figure saved as: "+outfile_annot)

def plot_LG_GFsize(orthogroups_dict, outfile_name = ""):
    """
    Plot the avg. gene family size on all linkage groups, with and without miniprot paralogs, and 
    showing the proportion of which are single-exon
    """

    all_LGs = [m for m in lg.LinkageGroup]
    cats = ["annot_intronfull", "annot_intronless", "mini_intronfull", "mini_intronless"]
    per_lg = {l: {c: [] for c in cats} for l in all_LGs}

    def count_LGid_GF_size(hog_LG_dict):
        size_LG_dict = {lg_ : 0 for lg_ in hog_LG_dict}
        for lg_ , tr_list in hog_LG_dict.items():
            size_LG_dict[lg_]+= len(tr_list)
        return size_LG_dict

    for HOG_id, orthogroup in orthogroups_dict.items():
        # print(f"\n-------{HOG_id}")
        counts = defaultdict(lambda: dict.fromkeys(cats, 0))

        for OG_member in orthogroup.members:

            ## check for multi-exon
            if OG_member.is_miniprot:
                multi = len(OG_member.gff_feature.child_ids_list) > 1
            else:
                multi = not OG_member.gff_feature.is_intronless()

            cat = ("mini" if OG_member.is_miniprot else "annot") + ("_intronfull" if multi else "_intronless")
            counts[OG_member.LG][cat] += 1
        
        for l, c in counts.items():
            for k in cats:
                per_lg[l][k].append(c[k])

        # if HOG_id== "N0.HOG0000004":
        #     break

    if outfile_name == "":
        return dict(per_lg)
    
    else:
        lw = 2
        fs = 25
        ymax_factor = 1.25
        plot_LGs = [l for l in all_LGs if l != lg.LinkageGroup.LGU]

        vals_dict = {c : [mean(per_lg[l][c]) for l in plot_LGs] for c in cats}
        colors = {
            "mini_intronfull" : "#676F54", # dusty olive
            "mini_intronless" : "#676F54", # dusty olive
            "annot_intronfull" : "#5F4B66", # vintage grape
            "annot_intronless" : "#5F4B66", # vintage grape
        }
        hatching = {
            "mini_intronfull" : "", 
            "mini_intronless" : "//",
            "annot_intronfull" : "", 
            "annot_intronless" : "//",
        }
        leg_labels = {
            "mini_intronfull" : "miniprot (multi-exon)",
            "mini_intronless" : "miniprot (single-exon)",
            "annot_intronfull" : "annotation (multi-exon)",
            "annot_intronless" : "annotation (single-exon)",
        }
        
        hatch_color = '#ffffff' # '#E2D4CA' #kind of eggshell white
        plt.rcParams['hatch.linewidth'] = lw  # default is 1.0
        plt.rcParams['hatch.color'] = hatch_color

        plt.rcParams['text.usetex'] = True # use \\textit{{{}}} for species names
        plt.rcParams['text.latex.preamble'] = r'\usepackage{sfmath} \renewcommand{\familydefault}{\sfdefault}'
        plt.rcParams['font.family'] = 'sans-serif'

        fig, ax = plt.subplots(1, 1, figsize=(20, 8))

        x_coord = range(len(plot_LGs))
        print(all_LGs)

        prev = [0 for x in x_coord]
        for category,vals in vals_dict.items():
            # print(f"{category} -- {vals}\n\n\n")
            ### TODO implement species division here when everything else works
            yval = [vals[i]/lg.LG_species_numbers[lg_] for i,lg_ in enumerate(plot_LGs)]
            
            yval_str =[f"{lg.LG_names[lg_]} - {yval[i]:.3f}" for i,lg_ in enumerate(plot_LGs)]
            print(f"--{category}-- \n{yval_str}\n")

            ax.bar(x_coord, yval, bottom=prev, label=leg_labels[category], color= colors[category], hatch=hatching[category])
            prev = [v+ prev[i] for i,v in enumerate(yval)]

        ax.tick_params(axis='x', labelsize=fs)
        ax.tick_params(axis='y', labelsize=fs)
        xtick_labels = [lg.LG_names[lg_].replace("Bruchini", "\\textit{{{Bruchini}}}") if "Bruchini" in lg.LG_names[lg_] else lg.LG_names[lg_].replace("Diorhabda", "\\textit{{{Diorhabda}}}") for lg_ in plot_LGs]
        ax.set_xticks(ticks = x_coord, labels = xtick_labels, rotation=90, fontsize=fs)
        ax.set_ylabel('mean GF size', fontsize=fs)

        ## legend
        plt.rcParams.update({'hatch.color': "#3f3832ff"})
        dashed_handle = mpatches.Patch(hatch = "//", alpha = 0.0)
        dashed_label = "proportion of genes \nthat are single-exon"
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1], labels[::-1], fontsize=fs, ncol=2, loc='upper left')
        ax.set_ylim(0,3)
    
        plt.tight_layout()
        plt.savefig(outfile_name, dpi = 300, transparent = True)# , bbox_inches='tight')
        print("Figure saved in the current working directory directory as: "+outfile_name)

        return dict(per_lg)

if __name__ == "__main__":
    
    username="milena"

    data_dir = f"/Users/{username}/work/PhD_code/PhD_chapter2/data/orthofinder"

    orthogroups_file = f"{data_dir}/N0.tsv" 
    unassigned_genes_path = f"{data_dir}/unassigned_genes.tsv" 
    annot_dict = data_paths.annotations_dict(username=username)
    miniprot_dict = data_paths.get_miniprot_paths(username=username)


    if True:
        species_order = [
            "A_obtectus",
            "C_maculatus",
            # "C_chinensis",
            "B_siliquastri",
            "B_varius",
            "D_carinulata",
            "D_sublineata"
        ]

        orthogroups = og.parse_orthogroups_class(
            filepath=orthogroups_file,
            annot_species=species_order,
            annotations_dict=annot_dict,
            add_gff_feature=True,
            unassigned_genes_path=unassigned_genes_path,
            miniprot_paths_dict=miniprot_dict)

        print(orthogroups["N0.HOG0000056"])
        
    ### plot Genome size vs. presence on linkage group scatter    
    if False:
        GS_vs_LG_number_correlations(orthogroups, species_order=species_order, filename_prefix=f"{data_dir}/translocations_sexchr_vs_A_numbers.png", max_GF_size=50)

    if True:
        # !!! remove duplicate cds from aobt bvar cchi cmac with 
        # sed '/^[^#\t]*\t[^\t]*\tCDS\t/d'
        
        c = plot_LG_GFsize(
            orthogroups_dict = orthogroups,
            outfile_name=f"{data_dir}/GFsize_on_linkagegroups.png"
        )
        # ## TODO no intronless annotations?
        ## TODO what exactly am i computing even, the mean GF size seems absurd???
        # print(c)
    if False:
        earliest_duplications, unassigned_geneIDs = og.parse_duplications(duplications_path=f"{data_dir}/Duplications.tsv", orthogroups_path=orthogroups_file, min_support=0.5)
