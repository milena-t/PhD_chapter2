"""
calculate gene density as bp_annotated_exon_or_cds / bp_total_from_fai
"""
from pairwise_wga_heatmap import fasta_indices
from sex_chromosomes import get_contig_names
from plot_gene_counts import get_miniprot_paths
import miniprot_stats_comparison as minialn
import parse_gff as gff
from statistics import mean
from scipy.stats import sem
import matplotlib.pyplot as plt
from circos.make_circos_karyotype_file import autosomes_lists
from matplotlib.ticker import FuncFormatter
import numpy as np
from GF_sizes import genome_sizes


def get_exon_coverage(miniprot_filepath, miniprot=True, gff_annot=False):
    """
    read the gff files and make dicts of {
        species : {contig : exon_bp_count, contig : exon_bp_count , ...}
    }
    """
    outdict = {}
    if gff_annot==True:
        annot_dict = gff.parse_gff3_general(filepath = miniprot_filepath, keep_feature_category=gff.FeatureCategory.Exon)
        # print(list(annot_dict.keys())[:50])
        miniprot=False
    if miniprot==True:
        annot_dict = minialn.miniprot_parse_alignment(miniprot_file=miniprot_filepath, include_cds=True)
    

    for miniID, gene_feature in annot_dict.items():
        try:
            outdict.setdefault(gene_feature.contig, 0) # if contig not yet in the dict then add it as new key with value 0
        except:
            print(f"\n{gene_feature}\n")
            raise RuntimeError
        outdict[gene_feature.contig] += gene_feature.length()

    return outdict



def calculate_gene_density(contig_lengths, species_gff, sex_chromosomes, autosomes, miniprot=False, get_gene_counts=False):

    if miniprot:
        # miniprot_exon_counts = get_exon_counts(species_miniprot, miniprot=True)
        gff_exon_counts = get_exon_coverage(species_gff, miniprot=True)
    else:
        gff_exon_counts = get_exon_coverage(species_gff, gff_annot=True)


    # gene_density = { contig : {"all" : 0.0, "gene" : 0.0, "mini" : 0.0} for contig in contig_lengths.keys() }
    gene_density = { contig : [] for contig in contig_lengths.keys() }
    if get_gene_counts:
        contig_annot = gff.parse_gff3_by_contig(species_gff)
        cds_features_counts = { contig : 0 for contig in contig_lengths.keys() }
        for contig, length in contig_lengths.items():
            try:
                cds_features_counts[contig] = len(contig_annot[contig])/length
            except:
                cds_features_counts[contig] = 0.0
        
    else:
        cds_features_counts = { contig : 0.0 for contig in contig_lengths.keys() }

    ## calculate densities for every contig
    for contig, contig_length in contig_lengths.items():
        try:
            gff_count = gff_exon_counts[contig]
        except:
            gff_count = 0
        
        gene_density[contig] = gff_count/contig_length
    
    ## calculate mean density for sex chromosomes, autosomes, and unplaced scaffolds
    sex_chromosomes["A"] = autosomes
    placed_scaffolds = sex_chromosomes["A"] + sex_chromosomes["X"] + sex_chromosomes["Y"]
    all_scaffolds = list(contig_lengths.keys())
    sex_chromosomes["unplaced"] = list(set(all_scaffolds) - set(placed_scaffolds))

    # densities_category = { category : {"all" : [], "gene" : [], "mini" : []} for category in sex_chromosomes.keys()}
    densities_category = { category : [] for category in sex_chromosomes.keys()}
    densities_count = { category : [] for category in sex_chromosomes.keys()}
    for chr_category, contig_list in sex_chromosomes.items():
        for contig in contig_list:
            densities_category[chr_category].append(gene_density[contig])
            densities_count[chr_category].append(cds_features_counts[contig])

    return densities_category, densities_count




def plot_gene_density(gene_density_dict, annot,  outfile_name = "", chromosome_categories = ["unplaced", "A", "X", "Y"], genome_sizes = {}):

    plt.rcParams['text.usetex'] = True # use \\textit{{{}}} for species names
    plt.rcParams['text.latex.preamble'] = r'\usepackage{sfmath} \renewcommand{\familydefault}{\sfdefault}'
    plt.rcParams['font.family'] = 'sans-serif'

    colors = {
        "A" : "#748B7A",
        "X" : "#BD351E", # red
        "Y" : "#5E79BD", # blue
        "unplaced" : "#909bb5"
    }
    legend_label = {
        "A" : "Autosomes",
        "unplaced" : "unpl. scaffolds",
        "X" : "X chromosome",
        "Y" : "Y chromosome"
    }
    
    fs = 30 # set font size
    species_names = sorted(list(gene_density_dict.keys()))
    
    point_offset=0.05
    xtick_ticklabel_pos = [i for i in range(len(species_names))]
    
    
    print(f" * {annot}")
    fig, ax = plt.subplots(1, 1, figsize=(13, 8))

    for i,cat in enumerate(chromosome_categories):
        mean_list = [0.0 for s in species_names]
        stderr_list = [0.0 for s in species_names]
        for j, species in enumerate(species_names):
            if gene_density_dict[species][cat]==[]:
                if species == "C_chinensis" and cat == "A":
                    mean_list[j] = mean(gene_density_dict[species]["unplaced"])
                    stderr_list[j] = sem(gene_density_dict[species]["unplaced"])
                else:
                    mean_list[j] = np.nan
                    stderr_list[j] = np.nan
            else:
                mean_list[j] = mean(gene_density_dict[species][cat])
                stderr_list[j] = sem(gene_density_dict[species][cat])
        xtick_pos = [k+point_offset*i for k in xtick_ticklabel_pos]
        # ax.plot(xtick_pos, mean_list, label = legend_label[cat], color = colors[cat], linewidth = 4) 
        print(f"\t{cat}: {mean_list}")
        ax.errorbar(xtick_pos, mean_list, yerr = stderr_list, color=colors[cat], linewidth =3, marker = ".", markersize=20, linestyle = ":", label = legend_label[cat])

    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: '' if x > 1 else f'{x*100.0:.0f}\%'))
    ylab = f"gene density"
    ax.set_ylabel(ylab, fontsize = fs)
    ax.tick_params(axis ='y', labelsize = fs)  
    # ymin,ymax = ax.get_ylim()
    # ax.set_ylim(0,ymax)

    ax.set_xticks([j+point_offset*1.5 for j in xtick_ticklabel_pos])
    species_axis_labels = [species.replace("_", ". ") for species in species_names]
    ax.set_xticklabels([f"\\textit{{{species}}}" for species in species_axis_labels], rotation=90, fontsize=fs)
    
    # set grid only for X axis ticks 
    ax.grid(True)
    ax.yaxis.grid(False)
    
    ax.tick_params(axis='y', labelsize=fs)
    ax.legend(fontsize=fs)

    plt.tight_layout()

    outfile_annot = outfile_name.replace(".png", f"_{annot}.png")
    plt.savefig(outfile_annot, dpi = 300, transparent = True)# , bbox_inches='tight')
    print("Figure saved as: "+outfile_annot)


def plot_gene_density_vs_GS(gene_density_dict, annot,  outfile_name = "", chromosome_categories = ["unplaced", "A", "X", "Y"], genome_sizes = {}, legend_loc="best"):

    ## get species names in order of size
    sizes_sorted = sorted(list(genome_sizes.values()))
    genome_sizes_ = {size : species for species,size in genome_sizes.items()} 
    species_names = [genome_sizes_[size] for size in sizes_sorted]

    plt.rcParams['text.usetex'] = True # use \\textit{{{}}} for species names
    plt.rcParams['text.latex.preamble'] = r'\usepackage{sfmath} \renewcommand{\familydefault}{\sfdefault}'
    plt.rcParams['font.family'] = 'sans-serif'

    colors_dict = {
        "A" : "#748B7A",
        "X" : "#BD351E", # red
        "Y" : "#5E79BD", # blue
        "Y_edge" : "#374C6E", # dusk blue darker
        "X_edge" : "#771C2C", # dark amaranth
        "A_edge" : "#376F59", # deep teal
        "unplaced" : "#909bb5",
        "vline" : '#9499A5', # cool steel
        "vtext" : "#5B6378", #blue slate
    }
    legend_label = {
        "A" : "Autosomes",
        "unplaced" : "unpl. scaffolds",
        "X" : "X chromosome",
        "Y" : "Y chromosome"
    }
    fs = 15 # font size
    ps = 20
    lw = 2
    
    point_offset=0.05
    xtick_ticklabel_pos = [i for i in range(len(species_names))]
    
    aspect_ratio = 18 / 14 # height / width
    height_pixels = 1400  # Height in pixels
    dpi = 300
    width_pixels = int(height_pixels * aspect_ratio)  # Width in pixels

    fig, ax = plt.subplots(1,1,figsize=(width_pixels/dpi, height_pixels/dpi))

    GS_list = [genome_sizes[s] for s in species_names]
    for x,s in zip(GS_list,species_names):
        ax.axvline(x=x, color=colors_dict["vline"], linestyle='-', linewidth=lw*0.25)
        species = s.replace("_", ". ")
        ax.text(
        x, 0.98, f"\\textit{{{species}}}",
        rotation=90,transform=ax.get_xaxis_transform(),  # x: data, y: axes fraction
        va="top",        # text hangs down from y=0.98
        ha="right",      # sits just left of the line (use "left" for right side)
        fontsize=fs*0.8, color=colors_dict["vtext"]
    )
    x_offset = {c:i*max(GS_list)/120 for c,i in zip(chromosome_categories,range(-1,2))}

    for i,cat in enumerate(chromosome_categories):
        gene_density = [0.0 for s in species_names]
        gd_errors = [0.0 for s in species_names]
        for j, species in enumerate(species_names):
            if gene_density_dict[species][cat]==[]:
                if species == "C_chinensis" and cat == "A":
                    gene_density[j] = mean(gene_density_dict[species]["unplaced"])
                    gd_errors[j] = sem(gene_density_dict[species]["unplaced"])
                else:
                    gene_density[j] = np.nan
                    gd_errors[j] = np.nan
            else:
                gene_density[j] = mean(gene_density_dict[species][cat])
                gd_errors[j] = sem(gene_density_dict[species][cat])
        x_coord = [i+x_offset[cat] for i in GS_list]
        # ax.plot(xtick_pos, mean_list, label = legend_label[cat], color = colors[cat], linewidth = 4) 
        print(f"\t{cat}: {gene_density}")

        ax.scatter(x_coord, gene_density, color = colors_dict[cat], s=ps, label = legend_label[cat])
        ax.errorbar(x_coord, gene_density, yerr=gd_errors, 
            # fmt="none", 
            color = colors_dict[cat], ecolor=colors_dict[f"{cat}_edge"], 
            capsize=4, elinewidth=lw, 
            linestyle = ":", linewidth =lw*0.5)
        # ax.errorbar(xtick_pos, mean_list, yerr = stderr_list, color=colors[cat], linewidth =3, marker = ".", markersize=20, linestyle = ":", label = legend_label[cat])

    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: '' if x > 1 else f'{x*100.0:.0f}\%'))
    ylab = f"gene density"
    ax.set_ylabel(ylab, fontsize = fs)
    ax.tick_params(axis ='y', labelsize = fs)  
    ax.tick_params(axis ='x', labelsize = fs)  
    ax.set_xlabel(f"Genome size in Mb", fontsize = fs)
    # ymin,ymax = ax.get_ylim()
    # ax.set_ylim(0,ymax)
    
    ax.tick_params(axis='y', labelsize=fs)
    ax.legend(fontsize=fs, loc=legend_loc)

    plt.tight_layout()

    outfile_annot = outfile_name.replace(".png", f"_{annot}.png")
    plt.savefig(outfile_annot, dpi = 300, transparent = True)# , bbox_inches='tight')
    print("Figure saved as: "+outfile_annot)


if __name__ == "__main__":
    username = "miltr339"
    faidx_dicts,_ = fasta_indices(username=username)
    miniprot_dict = get_miniprot_paths(username=username)
    annot_path_dict = minialn.get_gff_annot_paths(username=username)
    sex_chromosomes_dict = get_contig_names()
    autosomes_dict = autosomes_lists()
    data_dir = f"/Users/{username}/work/PhD_code/PhD_chapter2/data/"

    
    ## get gene densities for annotated genes from the gff
    gene_densities = {}
    gene_counts = {}

    plot_miniprot = True


    for species, miniprot_path in miniprot_dict.items():
        print(f"\n====================== {species} ======================")
        if plot_miniprot:
            annot_path = miniprot_path
            get_gene_counts = False
            annot="mini"
        else:
            annot_path = annot_path_dict[species]
            get_gene_counts = True
            annot="gene"
        
        gene_densities[species], gene_counts[species] = calculate_gene_density(
            contig_lengths=faidx_dicts[species], 
            species_gff=annot_path,
            sex_chromosomes=sex_chromosomes_dict[species],
            autosomes=autosomes_dict[species], 
            miniprot=plot_miniprot, get_gene_counts=get_gene_counts
            )
        
        if plot_miniprot:
            for chr_category, density_list in gene_densities[species].items():
                gene_count = gene_counts[species][chr_category]
                try:
                    print(f" - {chr_category} : {mean(density_list):.4f} bp annotated as exons")
                except:
                    print(f" - {chr_category} : NA, ({len(density_list)} contigs)")

        else:
            for chr_category, density_list in gene_densities[species].items():
                gene_count = gene_counts[species][chr_category]
                try:
                    print(f" - {chr_category} : {mean(density_list):.3f} bp annotated as exons ({len(density_list)} contigs), (avg. {mean(gene_count)*1000000:.3f} genes per Mb)")
                except:
                    print(f" - {chr_category} : NA, ({len(density_list)} contigs)")


    # plot_gene_density(
    #     gene_density_dict=gene_densities, annot=annot, 
    #     outfile_name=f"{data_dir}sex_chromosome_gene_density.png",
    #     chromosome_categories = ["A", "X", "Y"])

    if plot_miniprot:
        legend_loc = "center left"
    else:
        legend_loc = "center"

    plot_gene_density_vs_GS(
        gene_density_dict=gene_densities, annot=annot, # species_names=list(genome_sizes.keys()),
        outfile_name=f"{data_dir}sex_chromosome_gene_density_vs_GS.png",
        chromosome_categories = ["A", "X", "Y"],genome_sizes=genome_sizes, legend_loc=legend_loc)
