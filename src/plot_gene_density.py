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


def get_exon_counts(miniprot_filepath, miniprot=True, gff_annot=False):
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
        annot_dict = minialn.miniprot_parse_alignment(miniprot_file=miniprot_filepath, queryIDs=False)
    

    for miniID, gene_feature in annot_dict.items():
        outdict.setdefault(gene_feature.contig, 0) # if contig not yet in the dict then add it as new key with value 0
        outdict[gene_feature.contig] += gene_feature.length()

    return outdict


def calculate_gene_density(contig_lengths, species_miniprot, species_gff, sex_chromosomes, autosomes):

    miniprot_exon_counts = get_exon_counts(species_miniprot, miniprot=True)
    gff_exon_counts = get_exon_counts(species_gff, gff_annot=True)
    gene_density = { contig : {"all" : 0.0, "gene" : 0.0, "mini" : 0.0} for contig in contig_lengths.keys() }

    ## calculate densities for every contig
    for contig, contig_length in contig_lengths.items():
        try:
            gff_count = gff_exon_counts[contig]
        except:
            gff_count = 0
        try:
            mini_count = miniprot_exon_counts[contig]
        except:
            mini_count = 0
        
        all_count = gff_count+mini_count
        gene_density[contig]["all"] = all_count/contig_length
        gene_density[contig]["gene"] = gff_count/contig_length
        gene_density[contig]["mini"] = mini_count/contig_length
    
    ## calculate mean density for sex chromosomes, autosomes, and unplaced scaffolds
    sex_chromosomes["A"] = autosomes
    placed_scaffolds = sex_chromosomes["A"] + sex_chromosomes["X"] + sex_chromosomes["Y"]
    all_scaffolds = list(contig_lengths.keys())
    sex_chromosomes["unplaced"] = list(set(all_scaffolds) - set(placed_scaffolds))

    densities_category = { category : {"all" : [], "gene" : [], "mini" : []} for category in sex_chromosomes.keys()}
    for chr_category, contig_list in sex_chromosomes.items():
        for gene_category in ["all","gene","mini"]:
            for contig in contig_list:
                densities_category[chr_category][gene_category].append(gene_density[contig][gene_category])

    return densities_category

def plot_gene_density(gene_density_dict, outfile_name = ""):
    chromosome_categories = ["unplaced", "A", "X", "Y"]

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
    
    for annot in ["all","gene","mini"]:
        print(f" * {annot}")
        fig, ax = plt.subplots(1, 1, figsize=(13, 8))

        for i,cat in enumerate(chromosome_categories):
            mean_list = [0.0 for s in species_names]
            stderr_list = [0.0 for s in species_names]
            for j, species in enumerate(species_names):
                if gene_density_dict[species][cat][annot]==[]:
                    mean_list[j] = np.nan
                    stderr_list[j] = np.nan
                else:
                    mean_list[j] = mean(gene_density_dict[species][cat][annot])
                    stderr_list[j] = sem(gene_density_dict[species][cat][annot])
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


if __name__ == "__main__":
    username = "miltr339"
    faidx_dicts,_ = fasta_indices(username=username)
    miniprot_dict = get_miniprot_paths(username=username)
    annot_path_dict = minialn.get_gff_annot_paths(username=username)
    sex_chromosomes_dict = get_contig_names()
    autosomes_dict = autosomes_lists()
    data_dir = f"/Users/{username}/work/PhD_code/PhD_chapter2/data/"

    gene_densities = {}
    for species, miniprot_path in miniprot_dict.items():
        print(f"\n====================== {species} ======================")
        
        gene_densities[species] = calculate_gene_density(
            contig_lengths=faidx_dicts[species], 
            species_miniprot=miniprot_path,
            species_gff=annot_path_dict[species],
            sex_chromosomes=sex_chromosomes_dict[species],
            autosomes=autosomes_dict[species]
            )
        
        # for chr_category, density_dict in gene_densities[species].items():
        #     print(f" - {chr_category} : {density_dict}\n\n")


    plot_gene_density(gene_density_dict=gene_densities, outfile_name=f"{data_dir}sex_chromosome_gene_density.png")
