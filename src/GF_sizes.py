# check the position of the significant orthogroups in each annotation:
#  - do they look like tandem duplications or are they dispersed?
#  - Are the n basepairs around the gene enriched for a certain kind of TE compared to the rest of the genome?
#    (pick n like 1e5? pretty short probably)
#  - Are they intron-less?

import upsetplot
import numpy as np
import pandas as pd
from statistics import mean,median
from scipy.stats import sem
import matplotlib.pyplot as plt
import warnings
from typing import Literal

import parse_gff as gff
from sex_chromosomes import get_contig_names
import miniprot_stats_comparison as minialn
import plot_gene_counts as data_paths
from parse_orthogroups import  parse_orthogroups_dict


def parse_orthogroups_with_gff_class(filepath, annotations_dict, sex_chr_dict, unassigned_genes_path = "", miniprot_paths_dict = {}):
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
    annot_placed_dict = parse_orthogroups_dict(filepath=filepath, species_list=annot_species, OG_header="HOG")
    if unassigned_genes_path == "":
        annot_HOG_dict = annot_placed_dict
    else:
        annot_unplaced_dict = parse_orthogroups_dict(filepath=unassigned_genes_path, species_list=annot_species, OG_header="Orthogroup")
        annot_HOG_dict = annot_placed_dict | annot_unplaced_dict
        # print(annot_HOG_dict["OG0014394"])

    print(f"----- read gff annotations -----")
    annot_gff_dict = { species : gff.parse_gff3_general(annotations_dict[species], keep_feature_category=gff.FeatureCategory.Transcript, verbose=False) for species in annot_species}
    print(f"--------------------------------")

    hog_sexchr_dict = {}

    if miniprot_paths_dict != {}:
        print(f"--- read miniprot annotations ---")
        miniprot_dict = { species : minialn.miniprot_parse_alignment(miniprot_paths_dict[species]) for species in annot_species}
        print(f"---------------------------------")
        mini_X = {species : 0 for species in annot_species}
        mini_Y = {species : 0 for species in annot_species}

        singletons_with_miniprot_paralogs = 0
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
                        
                        if geneID in miniprot_dict[species]:
                            for miniprot_ID in miniprot_dict[species][geneID]:
                                mini_contig = miniprot_ID.contig
                            if mini_contig in sex_chr_dict[species]["X"]:
                                sex_chr_counts["X"] += 1
                                mini_X[species] += 1
                            elif mini_contig in sex_chr_dict[species]["Y"]:
                                sex_chr_counts["Y"] += 1
                                mini_Y[species] += 1
                            else:
                                sex_chr_counts["O"] += 1
                    
                    if HOG_id[:2] == "OG" and sum(list(sex_chr_counts.values()))>1:
                        singletons_with_miniprot_paralogs += 1
                        # print(f"{HOG_id}:{species} -> geneid {geneIDs}, {len(miniprot_dict[species][geneID])} miniprot paralog(s)")
                else:
                    # print(f"   - {species} : 0")
                    pass
                sexchr_counts_dict[species] = sex_chr_counts
            hog_sexchr_dict[HOG_id] = sexchr_counts_dict
        print(f"---<>---> mini-X paralogs : {mini_X}\n---<>---> mini-Y paralogs : {mini_Y} \n")
        print(f"---<>---> {singletons_with_miniprot_paralogs} singleton-OG miniprot paralogs")

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



def upset_HOG_sex_chromosomes(hog_sexchr_dict, chr_string = "Y", plot_filename = "upsetplot.png", min_intersection_size = 0):

    hog_sexchr_unnested = []
    autosome_excl = 0
    chr_linked = 0

    for HOG_id, sexchr_species_dict in hog_sexchr_dict.items():
        species_incl = []
        for species, sexchr_dict in sexchr_species_dict.items():
            allchr_present = 0
            for chr in chr_string:
                if sexchr_dict[chr] >0:
                    allchr_present += 1
            if len(chr_string) == allchr_present:
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

    warnings.filterwarnings("ignore", category=FutureWarning)
    warnings.filterwarnings("ignore", category=UserWarning)
    data = upsetplot.from_memberships(hog_sexchr_unnested)
    upsetplot.UpSet(data, subset_size="count", sort_by="cardinality", sort_categories_by="input", show_counts=True, min_subset_size=min_intersection_size).plot(fig=fig)

    if len(chr_string)==1:
        chr = chr_string
    elif len(chr_string)>1:
        chr = "X and Y"

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


def get_sexchr_hog_size(sexchr_dict):
    """
    Get the orthogroup size from a dictionary like this:
    HOG_ID : {
        species_1 : { "O" : n, "X" : n, "Y" : n },
        species_2 : { "O" : m, "X" : m, "Y" : m },
        ...
    }
    """
    size = 0
    for species, sexchr_numbers in sexchr_dict.items():
        for sexchr,num_paralogs in sexchr_numbers.items():
            size += num_paralogs
    return size


def get_sexchr_gf_size(sexchr_dict):
    """
    Get the gene family size in the species from a dictionary like this:
    species_2 : { "O" : m, "X" : m, "Y" : m }
    """
    size = 0
    for sexchr,num_paralogs in sexchr_dict.items():
        size += num_paralogs
    return size


def check_sex_linked_GF_size(hog_sexchr_dict, species_list, chr_string: Literal["X", "Y", "XY"], excl_chr:Literal["X", "Y",""] = "", outfile="boxplot.png", min_OG_size = 2, min_GF_size = 0, ymax_plot = 0, plot_type:Literal["medians_box", "means_bar"]="means_bar"):
    """
    Take sexchr size dict and make list of orthogroup sizes that are X-linked/Y-linked 
    if excl_chr is specified, then all orthogroups that have members on this chromosome are excluded from any analysis
    """
    chr_sizes = { c : {species : [] for species in species_list} for c in chr_string}
    # chr_alt_sizes = {species : [] for species in species_list}
    O_sizes = {species : [] for species in species_list}
    excl_OG = []

    if chr_string == "XY":
        plot_type = "means_bar"
    
    count_xy_linked = {s : [] for s in species_list}
    for HOG_id, sexchr_species_dict in hog_sexchr_dict.items():
        
        if get_sexchr_hog_size(sexchr_species_dict) <min_OG_size:
            excl_OG.append(HOG_id)
            continue
        
        for species, sexchr_dict in sexchr_species_dict.items():
            
            # if get_sexchr_gf_size(sexchr_dict=sexchr_dict) ==0:
            if get_sexchr_gf_size(sexchr_dict=sexchr_dict) <=min_GF_size:
                continue

            if excl_chr != "":
                if sexchr_dict[excl_chr] >0:
                    continue # skip other chromosome entirely, don't even add to O-sizes

            chr_count = {c : 0 for c in chr_string}
            A_count = {c : 0 for c in chr_string}
            for chr_ in chr_string:
                if sexchr_dict[chr_] > 0:
                    chr_count[chr_]+= 1
                elif sexchr_dict["O"] > 0: # only count size when the relevant orthogroup is larger than 0 in the current species
                    A_count[chr_]+=1
            
            if len(chr_string)>1:
                if len(chr_string) == sum(chr_count.values()):
                    count_xy_linked[species].append(get_sexchr_gf_size(sexchr_dict=sexchr_dict))
                    # print(f"{HOG_id}:{species}:{chr_count} : {sexchr_species_dict[species]}")
                    continue # gene on X and Y -> do not include
            
                if len(chr_string) == sum(A_count.values()):
                    # only added to A_count every time -> not in X and not in Y -> A exclusive
                    O_sizes[species].append(get_sexchr_gf_size(sexchr_dict=sexchr_dict))
            
                # if not A-exclusive and also not on both X and Y (see 'continue' above)...
                else:
                    for chr_,ccount in chr_count.items():
                        if ccount>0:
                            chr_sizes[chr_][species].append(get_sexchr_gf_size(sexchr_dict=sexchr_dict))

            else:
                if sum(chr_count.values())>0:
                    chr_sizes[chr_string][species].append(get_sexchr_gf_size(sexchr_dict=sexchr_dict))
                elif sum(A_count.values())>0:
                    O_sizes[species].append(get_sexchr_gf_size(sexchr_dict=sexchr_dict))

    print(f"  orthogroups per species that are excluded because they are X and Y linked at the same time:")
    for s, xy_count in count_xy_linked.items():
        if len(xy_count)>1:
            sem_s = sem(xy_count)
        else:
            sem_s = np.nan
        print(f"   - {s} ({len(xy_count)}) -->\t mean {mean(xy_count):.3f} [SEM: {sem_s:.3f}])")

    hog_excl = [og for og in excl_OG if "HOG" in og]
    print(f"  {len(excl_OG)} gene families excluded since they have a size < {min_OG_size}\n  {len(excl_OG)-len(hog_excl)} true unplaced genes and {len(hog_excl)} HOGs (not sigletons, but might be part of Cmac_C and only one other species making them effectively singletons for this analysis) ")
    # for hog_id in hog_excl:
    #     print_dict_ = hog_sexchr_dict[hog_id]
    #     print(f"{hog_id} : {print_dict_}")


    ## make dicts into one interleaved ones for easier plotting
    data_dict = {}

    print(f"\n---\nplotted numbers: ")
    if len(chr_string)==1:
        for species in species_list:
            chr_size = chr_sizes[chr_string][species]
            data_dict[f"\\textit{{{species}}}\n{chr_string} ({len(chr_size)} GFs)"] = chr_size
            O_size = O_sizes[species]
            data_dict[f"\\textit{{{species}}}\nA ({len(O_size)} GFs)"] = O_size
            # print(f" * {species} : GF_size means ({chr_string}: {mean(chr_size):.3f}) and (A {mean(O_size):.3f})")
            if plot_type == "means_bar":
                try:
                    print(f" * {species} : GF_size means ({chr_string}: {mean(chr_size):.3f} (length {len(chr_size)}), A: {mean(O_size):.3f} (length {len(O_size)}))")
                except:
                    print(f" * {species} : GF_size cant calculate means! ({chr_string}: length {len(chr_size)}), A: (length {len(O_size)}))")
                    raise RuntimeError
            else:
                print(f" * {species} : GF_size medians ({chr_string}: {median(chr_size):.3f}, A: {median(O_size):.3f})")
    elif len(chr_string)>1:
        for i, species in enumerate(species_list):
            for chr_ in chr_string:
                chr_size = chr_sizes[chr_][species]
                # data_dict[f"\\textit{{{species}}}\n{chr_} ({len(chr_size)} GFs)"] = chr_size
                data_dict[f"({len(chr_size)}) {chr_}:{i}"] = chr_size

            # print(f" * {species} : GF_size means ({chr_string}: {mean(chr_size):.3f}) and (A {mean(O_size):.3f})")
            try:
                print(f"  * {species}")
                for c in chr_string:
                    print(f"\t({c}: mean {mean(chr_size):.3f} ; length {len(chr_size)})")
            except:
                print(f"{species}")
                for c in chr_string:
                    print(f"\t({c}: length {len(chr_size)})")
                # print(f"\t(A: length {len(O_size)})")
                raise RuntimeError

            O_size = O_sizes[species]
            print(f"\t(A: mean {mean(O_size):.3f} ; length {len(O_size)})")
            #data_dict[f"\\textit{{{species}}}\nA ({len(O_size)} GFs)"] = O_size
            data_dict[f"({len(O_size)}) A:{i}"] = O_size

    ### make boxplot to show median size
    
    plt.rcParams['text.usetex'] = True # use \textit{} for species names
    plt.rcParams['text.latex.preamble'] = r'\usepackage{sfmath} \renewcommand{\familydefault}{\sfdefault}'
    plt.rcParams['font.family'] = 'sans-serif'
    # plt.rcParams['font.size'] = 16
    fs = 15 # font size
    lw = 2 # line width

    colors_dict = {
        "Y_fill" : "#495E83", # dusk blue
        "Y_edge" : "#374C6E", # dusk blue darker
        "Y_medians" : "#A7CCED", # icy blue
        "X_fill" : "#AB354A", # cherry rose
        "X_edge" : "#771C2C", # dark amaranth
        "X_medians" : "#EA9AA9", # cotton candy
        "O_fill" : "#51997C", # seagrass
        "O_edge" : "#376F59", # deep teal
        "O_medians" : "#8FCEB5", # perl_aqua
    }

    aspect_ratio = 20 / 12 # height / width
    height_pixels = 1400  # Height in pixels
    dpi = 300
    width_pixels = int(height_pixels * aspect_ratio)  # Width in pixels

    fig, ax = plt.subplots(1,1,figsize=(width_pixels/dpi, height_pixels/dpi))

    # tick_labels = [f"{key}\n({len(lists)} GFs)" for key,lists in data_dict.items()]
    tick_labels = [k.split(":")[0] for k in data_dict.keys()]
    lists = [means_list for means_list in data_dict.values()]

    if plot_type == "medians_box":

        width = 0.7
        pos_adjust=width*0.125

        tick_pos = [i+pos_adjust if i%2==1 else i-pos_adjust for i in range(1,len(tick_labels)+1)]
        bp = ax.boxplot(lists, positions=tick_pos, widths=width, patch_artist=True)   
        # set axis labels
        ax.set_ylabel("gene family size", fontsize=fs)

        ## modify boxplot colors
        if True:
            for i, box in enumerate(bp['boxes']):
                if i%2==0:
                    box.set(facecolor=colors_dict[f"{chr_string}_fill"], edgecolor=colors_dict[f"{chr_string}_edge"], linewidth=2)
                else:
                    box.set(facecolor=colors_dict["O_fill"], edgecolor=colors_dict["O_edge"], linewidth=2)
            for i, median_ in enumerate(bp['medians']):
                if i%2==0:
                    median_.set(color=colors_dict[f'{chr_string}_medians'], linewidth=lw)
                else:
                    median_.set(color=colors_dict['O_medians'], linewidth=lw)
            for i, whisker in enumerate(bp['whiskers']):
                # print(f"whisker: {i}")
                if i//2 % 2==0:
                    whisker.set(color=colors_dict[f'{chr_string}_edge'], linestyle='-',linewidth=lw)
                else:
                    whisker.set(color=colors_dict['O_edge'], linestyle='-',linewidth=lw)
            for i, cap in enumerate(bp['caps']):
                if i//2 % 2==0:
                    cap.set(color=colors_dict[f'{chr_string}_edge'],linewidth=lw)
                else:
                    cap.set(color=colors_dict['O_edge'],linewidth=lw)
            for i, flier in enumerate(bp['fliers']):
                if i%2==0:
                    flier.set(marker='.', markerfacecolor=colors_dict[f'{chr_string}_edge'], markeredgecolor=colors_dict[f'{chr_string}_edge'])
                else:
                    flier.set(marker='.', markerfacecolor=colors_dict['O_edge'], markeredgecolor=colors_dict['O_edge'])

    elif plot_type == "means_bar":

        if len(chr_string)>1:
            tick_pos = [i for i in range(1,len(tick_labels)+1)]
            box_width = 1/5
            pos_adjust=box_width*1.3
            for i in range(len(tick_labels)):
                if i%3 == 0:
                    tick_pos[i] =1+ i//3 -pos_adjust
                if i%3 == 1:
                    tick_pos[i] =1+ i//3 
                if i%3 == 2:
                    tick_pos[i] =1+ i//3 +pos_adjust
                print(f"{i} : {tick_pos[i]}")
        else:
            box_width = 0.7
            pos_adjust=box_width*0.125
            tick_pos = [i+pos_adjust if i%2==1 else i-pos_adjust for i in range(1,len(tick_labels)+1)]

        ymax_plot = 0
        colors_list = []
        errors_colors = []
        for data_key in tick_labels:
            if " A" in data_key:
                colors_list.append(colors_dict["O_fill"])
                errors_colors.append(colors_dict["O_edge"])
            else:
                for chr_ in chr_string:
                    if f" {chr_}" in data_key:
                        colors_list.append(colors_dict[f"{chr_}_fill"])
                        errors_colors.append(colors_dict[f"{chr_}_edge"])
        
        # set axis labels
        ax.set_ylabel("mean gene family size", fontsize=fs)
        means_lists = [mean(gf_sizes) for gf_sizes in lists]
        sem_lists = [sem(gf_sizes) for gf_sizes in lists]
        ax.bar(x=tick_pos, height=means_lists, yerr=sem_lists, width=box_width, color = colors_list)

        #custom errorbars for matching colors
        for xi, m, e, c in zip(tick_pos, means_lists, sem_lists, errors_colors):
            ax.errorbar(xi, m, yerr=e, fmt="none", ecolor=c, capsize=4, elinewidth=1.5)

        if len(chr_string)>1:
            fs_factor=0.8
            ax2 = ax.secondary_xaxis('bottom')
            ax2.set_xticks([i+1 for i in range(len(species_list))])
            species_list_ = [s.replace("_", ". ") for s in species_list]
            ax2.set_xticklabels([f"\\textit{{{s}}}" for s in species_list_], fontsize=fs*fs_factor, rotation=90)
            ax2.spines['bottom'].set_position(('outward', 70))   # 70 fo rbelow
            ax2.xaxis.set_ticks_position('none')
            ax2.spines['bottom'].set_visible(False)
            ax2.tick_params(axis='x', labelsize=fs*fs_factor)

    # ax.set_yscale('log')
    ax.set_xlabel("")
    ax.tick_params(axis='x', labelsize=fs) 
    ax.set_xticks(ticks = tick_pos, labels = tick_labels, fontsize=fs, rotation=90)
    ax.tick_params(axis='y', labelsize=fs)
    if chr_string == "XY":
        chr_string_ = "X or Y"
    else:
        chr_string_ = chr_string
    title = f"Sizes of A and {chr_string_}-linked gene families"
    if min_GF_size >0:
        title = f"{title} " + r"(min. GF size $\geq 2$)"
    ax.set_title(title, fontsize=fs)

    if ymax_plot>0:
        ax.set_ylim(0.5,ymax_plot)
    # layout rect=(left, bottom, right, top)

    ax.axhline(y=1, color='#B78F85', linestyle='--', linewidth=lw)
    plt.tight_layout()# rect=[0.0, 0.05, 1, 1])

    # transparent background
    plt.savefig(outfile, dpi = dpi, transparent = True)
    print(f"plot saved in current working directory as: {outfile}")




def orthogroups_summary_stats(hog_sexchr_dict = {}, filepath = "", annot_species = []):
    """
    calculate mean orthogroup and gene family size
    """
    if hog_sexchr_dict != {}:
        print(f"\n-------------------------------------------------------")
        print(f" results from orthogroups after sex chromosome assignment")
        test_og="OG0019090"
        # print(f"{test_og} : {hog_sexchr_dict[test_og]}")
        OG_sizes = {og : 0 for og in hog_sexchr_dict.keys()}
        GF_sizes = {s : [] for s in annot_species}
        for orthogroup, orthogroup_dict in hog_sexchr_dict.items():
            OG_size = 0
            for species, chr_dict in orthogroup_dict.items():
                gf_size = sum(list(chr_dict.values()))
                if gf_size == 0:
                    continue
                # if orthogroup[:2] == "OG" and gf_size>1:
                #     raise RuntimeError(f"parsing error for unplaced gene in orthogroup {orthogroup}: {orthogroup_dict}")
                GF_sizes[species].append(gf_size)
                OG_size += gf_size
            OG_sizes[orthogroup] = OG_size
        
        OG_sizes_list = list(OG_sizes.values())
        
        print(f"mean orthogroup size: {mean(OG_sizes_list):.3f} [SEM: {sem(OG_sizes_list):.3f}]")
        print(f"gene family sizes:")
        for species , gf_sizes in GF_sizes.items():
            print(f" * {species} ({len(gf_sizes)} gene families):\t mean: {mean(gf_sizes):.3f} [SEM: {sem(gf_sizes):.3f}]")

        print(f"\n")
    
    if filepath != "":
        assert len(annot_species)>0

        print(f"\n-----------------------------------------------------")
        filepath_ = filepath.split("/")[-1]
        print(f" results directly from {filepath_}, no sex chromosome assignment")
        
        orthogroups_dict = parse_orthogroups_dict(filepath=filepath, species_list=annot_species, OG_header="HOG")
        # print(orthogroups_dict["N0.HOG0000061"])
        OG_sizes = {og : 0 for og in orthogroups_dict.keys()}
        GF_sizes = {s : [] for s in annot_species}
        for orthogroup, orthogroup_dict in orthogroups_dict.items():
            OG_size = 0
            for species, list_ in orthogroup_dict.items():
                if isinstance(list_, str):
                    list_ = list_.split(", ")
                    gf_size = len(list_)
                    GF_sizes[species].append(gf_size)
                    OG_size += gf_size
                else:
                    list_ = []
                    gf_size = len(list_)
                OG_sizes[orthogroup] = OG_size

        OG_sizes_list = list(OG_sizes.values())
        print(f"mean orthogroup size: {mean(OG_sizes_list):.3f} [SEM: {sem(OG_sizes_list):.3f}]")
        print(f"gene family sizes:")
        for species , gf_sizes in GF_sizes.items():
            print(f" * {species} ({len(gf_sizes)} gene families):\t mean: {mean(gf_sizes):.3f} [SEM: {sem(gf_sizes):.3f}]")
                    





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

    #################################################
    only_gff = True    # if True: only read the orthogroups from orthofinder and don't add the within-species paralogs from miniprot ; if False: include the miniprot paralogs
    #################################################
    
    if only_gff:
        miniprot_dict = {}
        mini_filename = ""
    else:
        miniprot_dict = data_paths.get_miniprot_paths(username=username)
        mini_filename = "with_mini_paralogs_"

    if False:
        mini_hog_sexchr_dict = parse_orthogroups_with_gff_class(filepath=orthogroups_file, annotations_dict=annot_dict, unassigned_genes_path=unassigned_genes_path, sex_chr_dict=sex_chromosome_contigs, miniprot_paths_dict=miniprot_dict)
        orthogroups_summary_stats(hog_sexchr_dict=mini_hog_sexchr_dict, filepath=orthogroups_file, annot_species=species_order)
    

    ### plot the upsetplot for the different chromosome categories
    if False:
        ### read dict with sex chromosome 
        mini_hog_sexchr_dict = parse_orthogroups_with_gff_class(filepath=orthogroups_file, annotations_dict=annot_dict, unassigned_genes_path=unassigned_genes_path, sex_chr_dict=sex_chromosome_contigs, miniprot_paths_dict=miniprot_dict)
        
        min_intersection_size = {"X" : 20, "Y" : 0, "XY" : 0}
        for chr in ["Y","X","XY"]:
            print(f"\n\n>>>>> {chr} <<<<<")
            upset_HOG_sex_chromosomes(hog_sexchr_dict=mini_hog_sexchr_dict, chr_string=chr, plot_filename=f"{data_dir}/orthogroup_presence_{mini_filename}{chr}_upsetplot.png", min_intersection_size=min_intersection_size[chr])

    ### check if gene fmailies with X/Y members are on average larger
    if True:
        mini_hog_sexchr_dict = parse_orthogroups_with_gff_class(filepath=orthogroups_file, annotations_dict=annot_dict, unassigned_genes_path=unassigned_genes_path, sex_chr_dict=sex_chromosome_contigs, miniprot_paths_dict=miniprot_dict)
        orthogroups_summary_stats(hog_sexchr_dict=mini_hog_sexchr_dict, filepath=orthogroups_file, annot_species=species_order)

        if False:
            ymax_plot = {"X" : 20, "Y" : 0} # specify y limit in gene family size for the plot  (no filtering of the GF size data itself!)
            ## plot X and Y in separate plots
            for chr in ["Y","X"]:
                excl_chr = "X"
                if chr == "X":
                    excl_chr = "Y"
                filename = f"{data_dir}/GF_sizes_{chr}-linked_vs_A_comparison.png"

                print(f"\n\n>>>>> {chr} <<<<< (excl: {excl_chr})")
                check_sex_linked_GF_size(hog_sexchr_dict=mini_hog_sexchr_dict, species_list=species_order, chr_string=chr, excl_chr=excl_chr, outfile=filename, ymax_plot = ymax_plot[chr])
        if True:
            # plot three-color bar chart with all of them
            chr_ = "XY"
            if only_gff:
                filename = f"{data_dir}/GF_sizes_{chr_}-linked_vs_A_comparison_no_miniprot.png"
            else:
                filename = f"{data_dir}/GF_sizes_{chr_}-linked_vs_A_comparison.png"
            
            mingf = 2 # gene families need to have at least two members in the species
            if mingf>0:
                filename = filename.replace(".png", f"_gfsize_min{mingf}.png")

            print(f"\n\n>>>>> {chr_} <<<<< ")
            check_sex_linked_GF_size(hog_sexchr_dict=mini_hog_sexchr_dict, species_list=species_order, chr_string=chr_, min_GF_size = mingf, outfile=filename, plot_type="means_bar")