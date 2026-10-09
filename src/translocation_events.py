import matplotlib.pyplot as plt
from statistics import mean
from scipy.stats import sem
from collections import defaultdict

import linkage_groups as lg
import GF_sizes as OGs
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
                X_nums.setdefault(GF_size, [LG_number]).append(LG_number)
            elif lg.LinkageGroup.LGX not in linkage_groups and lg.LinkageGroup.LGY in linkage_groups:
                Y_nums.setdefault(GF_size, [LG_number]).append(LG_number)
            elif lg.LinkageGroup.LGX in linkage_groups and lg.LinkageGroup.LGY in linkage_groups:
                XY_nums.setdefault(GF_size, [LG_number]).append(LG_number)
            else:
                A_nums.setdefault(GF_size, [LG_number]).append(LG_number)
        
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

def plot_LG_GFsize(orthogroups_dict, outfile_name):
    """
    Plot the avg. gene family size on all linkage groups, with and without miniprot paralogs, and 
    showing the proportion of which are single-exon
    """

    mini_intronless = defaultdict(list)
    mini_intronfull = defaultdict(list)
    annot_intronless = defaultdict(list)
    annot_intronfull = defaultdict(list)

    def count_LGid_GF_size(hog_LG_dict):
        size_LG_dict = {lg_ : 0 for lg_ in hog_LG_dict.keys()}
        for lg_ , tr_list in hog_LG_dict.items():
            size_LG_dict[lg_]+= len(tr_list)
        return size_LG_dict

    for HOG_id, orthogroup in orthogroups_dict.items():
        # print(f"\n-------{HOG_id}")

        hog_mini_intronless=defaultdict(int)
        hog_mini_intronfull=defaultdict(int)
        hog_annot_intronless=defaultdict(int)
        hog_annot_intronfull=defaultdict(int)

        for OG_member in orthogroup.members:
            feature_kids=OG_member.gff_feature.child_ids_list
            assert feature_kids !=None
            if OG_member.is_miniprot:
                assert feature_kids !=None
                if len(feature_kids)>1:
                    hog_mini_intronfull.setdefault(OG_member.LG, [OG_member.transcript_ID]).append(OG_member.transcript_ID)
                else:
                    print(OG_member)
                    hog_mini_intronless.setdefault(OG_member.LG, [OG_member.transcript_ID]).append(OG_member.transcript_ID)

            else:
                if len(feature_kids)>1:
                    hog_annot_intronfull.setdefault(OG_member.LG, [OG_member.transcript_ID]).append(OG_member.transcript_ID)
                else:
                    hog_annot_intronless.setdefault(OG_member.LG, [OG_member.transcript_ID]).append(OG_member.transcript_ID)
            
            # print(f" - {OG_member.transcript_ID} child features: {len(feature_kids)}: {feature_kids[:30]}")
        hog_mini_intronless = count_LGid_GF_size(hog_mini_intronless)
        hog_mini_intronfull = count_LGid_GF_size(hog_mini_intronfull)
        hog_annot_intronless = count_LGid_GF_size(hog_annot_intronless)
        hog_annot_intronfull = count_LGid_GF_size(hog_annot_intronfull)

        for lg_, gfsize in hog_mini_intronless.items():
            mini_intronless.setdefault(lg_, []).append(gfsize)
        for lg_, gfsize in hog_mini_intronfull.items():
            mini_intronfull.setdefault(lg_, []).append(gfsize)
        for lg_, gfsize in hog_annot_intronless.items():
            annot_intronless.setdefault(lg_, []).append(gfsize)
        for lg_, gfsize in hog_annot_intronfull.items():
            annot_intronfull.setdefault(lg_, []).append(gfsize)

    
    
    return dict(mini_intronless),dict(mini_intronfull),dict(annot_intronless),dict(annot_intronfull)


if __name__ == "__main__":
    
    username="milena"

    data_dir = f"/Users/{username}/work/PhD_code/PhD_chapter2/data/orthofinder"

    orthogroups_file = f"{data_dir}/N0.tsv" 
    unassigned_genes_path = f"{data_dir}/unassigned_genes.tsv" 
    annot_dict = data_paths.annotations_dict(username=username)
    miniprot_dict = data_paths.get_miniprot_paths(username=username)
    species_order = [
        "A_obtectus",
        "C_maculatus",
        "C_chinensis",
        "B_siliquastri",
        "B_varius",
        "D_carinulata",
        "D_sublineata"
    ]

    if True:
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
        mini_intronless,mini_intronfull,annot_intronless,annot_intronfull = plot_LG_GFsize(orthogroups_dict = orthogroups, outfile_name=f"{data_dir}/GFsize_on_linkagegroups.png")
        print(
            f"""
mini_intronless : {mini_intronless.keys()}
mini_intronfull : {mini_intronfull.keys()}
annot_intronless : {annot_intronless.keys()}
annot_intronfull : {annot_intronfull.keys()}
"""
        )

    if False:
        earliest_duplications, unassigned_geneIDs = og.parse_duplications(duplications_path=f"{data_dir}/Duplications.tsv", orthogroups_path=orthogroups_file, min_support=0.5)
