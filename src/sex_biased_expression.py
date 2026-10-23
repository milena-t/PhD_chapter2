"""
check for median differences in log2FC between genes sorted by different categories of gene families they belong to
"""


import parse_orthogroups as og
import parse_gff as gff
from sex_chromosomes import get_contig_names
import plot_gene_counts as data_paths
from GF_sizes import parse_orthogroups_with_gff_class

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pandas as pd


def logFC_files_dict(username="miltr339"):
    datadir = f"/Users/{username}/work/PhD_code/PhD_chapter2/data/DE_analysis"
    outdict = {
        "ab" : f"{datadir}/Cmac_sex_DE_abdomen_edgeR.txt",
        "ht" : f"{datadir}/Cmac_sex_DE_head_thorax_edgeR.txt",
    }
    return outdict


def plot_logFC_change(hog_sexchr_dict, parsed_gff_dict, logFC_path, N0_file, species_sexchr_dict, sexchr = "X", outfile = "", species = "C_maculatus",title="", yl=3):
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
    df = pd.read_csv(logFC_path, sep="\t")
    logfc_vals = df["logFC"]

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
                    try:
                        lfcval = logfc_vals[gID]
                    except:
                        continue # if gene not expressed

                    if hog_category == "A" and tr_contig == "A":
                        gene_class_dict["A_in_A_gf"].append(lfcval)

                    if hog_category == f"{sexchr}_all" and tr_contig == "A":
                        gene_class_dict[f"A_in_{sexchr}_all_gf"].append(lfcval)
                    if hog_category == f"{sexchr}_all" and tr_contig == "X":
                        gene_class_dict[f"X_in_{sexchr}_all_gf"].append(lfcval) 
                    if hog_category == f"{sexchr}_{species}" and tr_contig == "A":
                        gene_class_dict[f"A_in_{sexchr}_{species}_gf"].append(lfcval)
                    if hog_category == f"{sexchr}_{species}" and tr_contig == "X":
                        gene_class_dict[f"X_in_{sexchr}_{species}_gf"].append(lfcval) 
                    if hog_category ==  f"{sexchr}{other_sexchr}_{species}" and tr_contig == "A":
                        gene_class_dict[f"A_in_{sexchr}{other_sexchr}_all_gf"].append(lfcval)
                    if hog_category ==  f"{sexchr}{other_sexchr}_{species}" and tr_contig == "X":
                        gene_class_dict[f"X_in_{sexchr}{other_sexchr}_{species}_gf"].append(lfcval) 

    for c,l in gene_class_dict.items():
        print(f"{c} : {len(l)} genes")
    print(f"---------------------------------")

    ### make boxplot to show median size
    if True:

        colors_dict = {
            "Y_fill" : "#495E83", # dusk blue
            "Y_edge" : "#374C6E", # dusk blue darker
            "Y_medians" : "#A7CCED", # icy blue
            "X_fill" : "#AB354A", # cherry rose
            "X_edge" : "#771C2C", # dark amaranth
            "X_medians" : "#EA9AA9", # cotton candy
            "O_fill" : "#77667A", # dusty lavender
            "O_edge" : "#574A59", # mauve shadow
            "O_medians" : "#9E90A2", # lilac ash
        }
        
        species_ = species.replace("_", ". ")
        ticklabels_dict = {
            "A_in_A_gf" : f"exclusively\nAutosomal",
            
            f"A_in_{sexchr}_all_gf" : f"A-genes in conserved\n{sexchr}-linked GFs",
            f"X_in_{sexchr}_all_gf" : f"X-genes in conserved\n{sexchr}-linked GFs",
            
            f"A_in_{sexchr}_{species}_gf" : f"A-genes in \n\\textit{{{species_}}}-only\n{sexchr}-linked GFs",
            f"X_in_{sexchr}_{species}_gf" : f"X-genes in \n\\textit{{{species_}}}-only\n{sexchr}-linked GFs",

            f"A_in_{sexchr}{other_sexchr}_all_gf" : f"A-genes in {sexchr}\nand {other_sexchr}-linked GFs\nin all species",
            f"X_in_{sexchr}{other_sexchr}_{species}_gf" : f"X-genes in {sexchr}\nand {other_sexchr}-linked GFs in\n\\textit{{{species_}}}",
        }

        colors_data_assign = {
            "A_in_A_gf" : "O",
            f"A_in_{sexchr}_all_gf" : f"O",
            f"X_in_{sexchr}_all_gf" : sexchr,
            f"A_in_{sexchr}_{species}_gf" : f"O",
            f"X_in_{sexchr}_{species}_gf" : sexchr,
            f"A_in_{sexchr}{other_sexchr}_all_gf" : f"O",
            f"X_in_{sexchr}{other_sexchr}_{species}_gf" : sexchr,
        }

        plt.rcParams['text.usetex'] = True # use \textit{} for species names
        plt.rcParams['text.latex.preamble'] = r'\usepackage{sfmath} \renewcommand{\familydefault}{\sfdefault}'
        plt.rcParams['font.family'] = 'sans-serif'
        # plt.rcParams['font.size'] = 16
        fs = 15 # font size
        lw = 2 # line width

        aspect_ratio = 20 / 14 
        height_pixels = 1400  # Height in pixels
        dpi = 300
        width_pixels = int(height_pixels * aspect_ratio)  # Width in pixels

        fig, ax = plt.subplots(1,1,figsize=(width_pixels/dpi, height_pixels/dpi))
        # ax.yaxis.set_major_locator(MaxNLocator(integer=True))

        width = 0.7
        pos_adjust=width*0.125
        
        incl_samples = [f"{l}" for l,v in gene_class_dict.items() if len(v)>0]
        # tick_pos = [i+pos_adjust if i%2==1 else i-pos_adjust for i in range(1,len(xticklabels)+1)]
        tick_pos = range(1,len(incl_samples)+1)
        lists = [gene_class_dict[xtick] for xtick in incl_samples]
        xticklabels = [f"{ticklabels_dict[xtick]} ({len(gene_class_dict[xtick])})" for xtick in incl_samples if len(gene_class_dict[xtick])>0]

        bp = ax.boxplot(lists, positions=tick_pos, widths=width, patch_artist=True)   
        # set axis labels
        ax.set_ylabel("log2FC (F-M)", fontsize=fs)
        # ax.set_xlabel("categories", fontsize = fs, labelpad = 25)
        ax.tick_params(axis='x', labelsize=fs) 
        ax.set_xticks(ticks = tick_pos, labels = xticklabels, fontsize=fs, rotation=90)
        ax.tick_params(axis='y', labelsize=fs)
        ax.axhline(y=0, color='#B78F85', linestyle='--', linewidth=lw*0.75)
        ax.set_ylim(yl*-1,yl)

        ## modify boxplot colors
        if True:
            for i, box in enumerate(bp['boxes']): 
                sexchr = colors_data_assign[incl_samples[i]]   
                box.set(facecolor=colors_dict[f"{sexchr}_fill"], edgecolor=colors_dict[f"{sexchr}_edge"], linewidth=2)
            for i, median_ in enumerate(bp['medians']):
                sexchr = colors_data_assign[incl_samples[i]]
                median_.set(color=colors_dict[f'{sexchr}_medians'], linewidth=lw)
            for j, whisker in enumerate(bp['whiskers']):
                # print(f"whisker: {i}")
                i= j//2
                sexchr = colors_data_assign[incl_samples[i]]
                whisker.set(color=colors_dict[f'{sexchr}_edge'], linestyle='-',linewidth=lw)
            for j, cap in enumerate(bp['caps']):
                i= j//2
                sexchr = colors_data_assign[incl_samples[i]]
                cap.set(color=colors_dict[f'{sexchr}_edge'],linewidth=lw)
            for i, flier in enumerate(bp['fliers']):
                sexchr = colors_data_assign[incl_samples[i]]
                flier.set(marker='.', markerfacecolor=colors_dict[f'{sexchr}_edge'], markeredgecolor=colors_dict[f'{sexchr}_edge'])
        
        if title!="":
            plt.title(title, fontsize=fs*1.25)
        plt.tight_layout()# rect=[0.0, 0.05, 1, 1])
        # transparent background
        plt.savefig(outfile, dpi = dpi, transparent = True)
        print(f"plot saved in current working directory as: {outfile}")


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
        focal_species_ = focal_species.replace("_", ". ")
        tissue_fullname = {
            "ht" : "somatic",
            "ab" : "reproductive"
        }
        ylimit = {
            "ht" : 6,
            "ab" : 6
        }
        
        species_full_GFFdict = gff.parse_gff3_general(filepath=annot_dict[focal_species], verbose=False)

        for tissue in ["ht", "ab"]:
            for sexchr in ["X", "Y"]:
                print(f"\n\n>>>>>>>>>>>>>>> {sexchr} : {tissue_fullname[tissue]} <<<<<<<<<<<<<<<")
                plot_logFC_change(
                    hog_sexchr_dict=hog_sexchr_dict, 
                    parsed_gff_dict=species_full_GFFdict, 
                    logFC_path=logfc_files[tissue],
                    species_sexchr_dict = sex_chromosome_contigs[focal_species],
                    sexchr = sexchr, species=focal_species,
                    N0_file = orthogroups_file,
                    outfile=f"{outdir}/logFC_boxplot_{sexchr}_{tissue}.png",
                    title = f"Average expression sex-bias\nin \\textit{{{focal_species_}}} {tissue_fullname[tissue]} tissue",
                    yl=ylimit[tissue]
                )
        