"""
Calculate gene density (via exon annotated regions)
Also single exon gene proportions to compare miniprot and native annotations
"""

import sex_chromosomes
import parse_gff as gff
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

def plot_single_exon_no_species_specific_three_annot(native_numbers, miniprot_numbers, species_names, filename = ""):

    print(f" plotting for these {len(species_names)} species: \n{species_names}")

    # X coordinates for the groups
    x = np.arange(len(species_names))

    # figure proportions according to the data included (longer or shorter)

    # fontsize scales with the dpi somehow which i have to do extra because i change the aspect ratio manually below
    fs = 60 # 37 originally
    plt.rcParams['text.usetex'] = True # use \textit{} for species names
    plt.rcParams['text.latex.preamble'] = r'\usepackage{sfmath} \renewcommand{\familydefault}{\sfdefault}'
    plt.rcParams['font.family'] = 'sans-serif'

    width = 0.33 
    aspect_ratio = 21 / 12 
    ymax_factor = 1.25

    height_pixels = 2000  # Height in pixels
    width_pixels = int(height_pixels * aspect_ratio)  # Width in pixels
    fig = plt.figure(figsize=(width_pixels / 100, height_pixels / 100), dpi=100)
    ax = fig.add_subplot(111)


    colors = {
        "native" : "#44506E", 
        "minimap" : "#8390B1",
    }
    hatch_color = '#ffffff' # '#E2D4CA' #kind of eggshell white
    plt.rcParams['hatch.color'] = hatch_color
    plt.rcParams['hatch.linewidth'] = 2  # default is 1.0
    hatch_label = f"proportion of which are intronless"

    #### plot native annotation ####
    num_single_exon_genes = native_numbers["num_single_exon_genes"]
    num_total_genes = native_numbers["num_total_genes"]

    single_exon_genes = [num_single_exon_genes[species] for species in species_names]
    multi_exon_genes = [num_total_genes[species]-num_single_exon_genes[species] for species in species_names]

    color = [colors["native"], colors["native"]]
    hatching = ["//", "" ]

    category = " (genome annotation)"
    x_subtr = width/2
    
    print(f"bar width = {width}")

    # total number of genes from annotations (with single-exons hatched)
    native_rects1_base = ax.bar(x - x_subtr, single_exon_genes, width, label=hatch_label, color= color[0], hatch=hatching[0])
    native_rects1_top = ax.bar(x - x_subtr, multi_exon_genes, width, bottom=single_exon_genes, label='all genes'+category, color= color[1], hatch=hatching[1])

    #### plot miniprot gene structures ####
    num_single_exon_genes = miniprot_numbers["num_single_exon_genes"]
    num_total_genes = miniprot_numbers["num_total_genes"]
    single_exon_genes = [num_single_exon_genes[species] for species in species_names]
    multi_exon_genes = [num_total_genes[species]-num_single_exon_genes[species] for species in species_names]

    color = [colors["minimap"], colors["minimap"]]
    hatching = ["//", ""]

    # total number of genes (with single-exons hatched)            
    orthodb_rects1_base = ax.bar(x + x_subtr, single_exon_genes, width, label=hatch_label, color= color[0], hatch=hatching[0])            
    orthodb_rects1_top = ax.bar(x + x_subtr, multi_exon_genes, width, bottom=single_exon_genes, label='miniprot gene structures', color= color[1], hatch=hatching[1])  

    plt.rcParams.update({'hatch.color': hatch_color})
    ymax_factor = 1.3
    ymax = max(num_total_genes.values())*ymax_factor

    #### set up labels and stuff ####
    
    ax.set_ylabel('Number of genes', fontsize=fs+4)
    ax.set_title('Proportion of single-exon genes', fontsize=fs+4)
    ax.set_xticks(x)
    xtick_labels = [species.replace("_", ". ") for species in species_names]
    ax.set_xlabel('', fontsize=fs+4)
    ax.set_xticklabels([f"\\textit{{{species}}}" for species in xtick_labels], rotation=90, fontsize=fs)
    ax.set_yticklabels([f'{int(tick)/1e3:.0f}k' for tick in ax.get_yticks()], fontsize=fs)

    # make custom legend patch for the dashed bars
    plt.rcParams.update({'hatch.color': "#3f3832ff"})
    dashed_handle = mpatches.Patch(hatch = "//", alpha = 0.0)
    dashed_label = "proportion of genes \nthat are single-exon"

    # Legend with custom order
    handles, labels = ax.get_legend_handles_labels()

    new_order = [1,3]
    handles = [handles[idx] for idx in new_order]
    labels = [labels[idx] for idx in new_order]
    handles.append(dashed_handle)
    labels.append(dashed_label)

    ax.legend(handles, labels, fontsize=fs, ncol=2, loc='upper center')

    # add space at the top of the plot for the legend
    ax.set_ylim(0, int(ymax))
    ax.set_xlim(-0.5, len(xtick_labels)-0.5)

    plt.tight_layout()

    plt.savefig(filename, dpi = 300, transparent = True)
    print("Figure saved in the current working directory directory as: "+filename)



def get_single_exon_genes(paths_dict, write_to_file = False, outfile_name = "", include_total_gene_num = False):
    """
    get a dictionary with {species_name : [list, of, single, exon, transcripts]} from files generated in 04c_make_single_exon_gene_lists.sh
    """
    species_dict = {}
    num_single_exon_genes = {}
    num_transcripts = {}

    all_break = False
    for species,infile_path in paths_dict.items():

        print_statements = False
        infile_dict = {}
        list_key = ""
        num_key = ""
        transcripts_num_key = ""
        with open(infile_path, "r") as infile:
            dict_elements = infile.readlines()
            for species_line in dict_elements:
                try:
                    key = species_line.split(":")[0]
                    value = ":".join(species_line.split(":")[1:]) # some of the transcript names have ":" in them to spite me personally
                except:
                    print_statements = True
                    print(species_line[0:150])
                    continue
                value_list = [gene.split("|")[-1] for gene in value.strip().split(",")] # TODO # native annotations have some stuff going on before "|" and to actually make it work you need to get rid of that
                infile_dict[key]=value_list
                if "list" in key:
                    list_key = key
                if "single exon transcripts" in key:
                    num_key = key
                if "number of transcripts" in key:
                # if "gene features" in key:    
                    transcripts_num_key = key
                
 
            if all_break:
                break
            if print_statements:
                print(species)
                print(list_key)
                print(infile_dict.keys())
                print_statements = False
        # species_dict[species] = infile_dict["list of single-exon transcript IDs"]
        try:
            species_dict[species] = infile_dict[list_key]
        except:
            print(f"{species} single-exon list did not work with {list_key}")
        
        try:
            num_single_exon_genes[species] = int(infile_dict[num_key][0])
        except:
            print(f"{species} single-exon number did not work with {num_key}")

        if include_total_gene_num:
            try:
                num_transcripts[species] = int(infile_dict[transcripts_num_key][0])
            except:
                keys = list(infile_dict.keys())
                print(f"{species} transcript number did not work with {transcripts_num_key}, available keys are: {keys}")

    # num_single_exon_genes = {species : len(IDs) for species, IDs in species_dict.items()}
    if not include_total_gene_num:
        return(species_dict, num_single_exon_genes)
    elif include_total_gene_num:
        return(species_dict, num_single_exon_genes, num_transcripts)


def get_gff_annot_paths(username = "miltr339"):
    filedir = f"/Users/{username}/work/chapter2/native_annotations/"
    outdict = {
        "A_obtectus" : f"{filedir}A_obtectus.gff",
        "B_siliquastri" : f"{filedir}B_siliquastri.gff",
        "B_varius" : f"{filedir}B_varius.gff",
        "C_chinensis" : f"{filedir}C_chinensis.gff",
        "C_maculatus" : f"{filedir}C_maculatus.gff",
        "D_carinulata" : f"{filedir}D_carinulata.gff",
        "D_sublineata" : f"{filedir}D_sublineata.gff",
    }
    return outdict

def get_miniprot_annot_paths(username="miltr339"):
    filedir = f"/Users/{username}/work/chapter2/miniprot_annot/"
    outdict = {
        "A_obtectus" : f"{filedir}A_obtectus_miniprot_no_cross_hits.gff",
        "B_siliquastri" : f"{filedir}B_siliquastri_miniprot_no_cross_hits.gff",
        "B_varius" : f"{filedir}B_varius_miniprot_no_cross_hits.gff",
        "C_chinensis" : f"{filedir}C_chinensis_miniprot_no_cross_hits.gff",
        "C_maculatus" : f"{filedir}C_maculatus_miniprot_no_cross_hits.gff",
        "D_carinulata" : f"{filedir}D_carinulata_miniprot_no_cross_hits.gff",
        "D_sublineata" : f"{filedir}D_sublineata_miniprot_no_cross_hits.gff",
    }
    return outdict

def single_exon_paths(username="miltr339"):
    filedir = f"/Users/{username}/work/PhD_code/PhD_chapter2/data/single_exon_stats"
    outdict_native = {
        "A_obtectus" : f"{filedir}/A_obtectus_single_exon_stats_with_transcript_list.txt",
        "B_siliquastri" : f"{filedir}/B_siliquastri_single_exon_stats_with_transcript_list.txt",
        "B_varius" : f"{filedir}/B_varius_single_exon_stats_with_transcript_list.txt",
        "C_chinensis" : f"{filedir}/C_chinensis_single_exon_stats_with_transcript_list.txt",
        "C_maculatus" : f"{filedir}/C_maculatus_single_exon_stats_with_transcript_list.txt",
        "D_carinulata" : f"{filedir}/D_carinulata_single_exon_stats_with_transcript_list.txt",
        "D_sublineata" : f"{filedir}/D_sublineata_single_exon_stats_with_transcript_list.txt",
    }
    outdict_miniprot = {
        "A_obtectus" : f"{filedir}/A_obtectus_miniprot_single_exon_stats_with_transcript_list.txt",
        "B_siliquastri" : f"{filedir}/B_siliquastri_miniprot_single_exon_stats_with_transcript_list.txt",
        "B_varius" : f"{filedir}/B_varius_miniprot_single_exon_stats_with_transcript_list.txt",
        "C_chinensis" : f"{filedir}/C_chinensis_miniprot_single_exon_stats_with_transcript_list.txt",
        "C_maculatus" : f"{filedir}/C_maculatus_miniprot_single_exon_stats_with_transcript_list.txt",
        "D_carinulata" : f"{filedir}/D_carinulata_miniprot_single_exon_stats_with_transcript_list.txt",
        "D_sublineata" : f"{filedir}/D_sublineata_miniprot_single_exon_stats_with_transcript_list.txt",
    }
    return outdict_native,outdict_miniprot


class MiniAln:
    """
    read miniprot alignment data 
    """
    def __init__(self, ID:str, target:str, rank:int, identity:float, contig:str, start:int, end:int) -> None:
        self.ID=ID
        self.target=target
        self.rank=rank
        self.start=start
        self.end=end
        self.contig=contig
        if identity<=1: # use percent not proportion
            self.identity=100*identity
        else:
            self.identity=identity
    def __repr__(self):
        return "MiniAln"
    def __str__(self) -> str:
        return f"""Miniprot Alignment ID: {self.ID}
 * query ID: {self.target}
 * aligned on contig: {self.contig} with {self.identity}% sequence identity
 * rank: {self.rank}"""
    def length(self):
        return abs(self.start-self.end)




def miniprot_parse_alignment(miniprot_file, queryIDs = True):
    """
    parse the miniprot alignment into a dictionary by query transcript ID
    {   
        target (genomeAnnot transcriptID) : [ MiniAln(class)1, MiniAln(class)2, ... ],
    }
    if queryIDs=False, then the dict is not nested with a list, and its jus the miniprot aln IDs with the class as the key
    """
    geneIDs_map_counts = {}
    count_paf = 0
    count_dup_paf =0
    count_aln = 0
    with open(miniprot_file, "r") as miniprot_infile:
        for mini_line in miniprot_infile.readlines(): # skip gff3 header
            if mini_line[0]=="#":
                continue
                line = mini_line.strip().split()
                count_paf+=1
                protein_target = line[1]
                if protein_target not in geneIDs_map_counts:
                    geneIDs_map_counts[protein_target] = [] 
                else:
                    count_dup_paf +=1
            else:
                line = mini_line.strip().split("\t")
                try:
                    contig,source,category,start,stop,score,strandedness,frame,attributes_=[c for c in line if len(c)>0]
                except Exception as e:
                    print(f"mini line could not be parsed! line: \n{mini_line}\nlist ({len(line)} items, should be 9)\n{line}\nerror:\n{e}")
                if category != "mRNA":
                    continue
                attributes={}
                for attr in attributes_.strip().split(";"):
                    key,value=attr.strip().split("=")[-2:]
                    if " " in value:
                        attributes[key]=value.split()[0]
                    else:
                        attributes[key]=value
                int(start),int(stop)
                mini_alignment = MiniAln(
                    ID=attributes["ID"], target=attributes["Target"], rank=int(attributes["Rank"]), identity=float(attributes["Identity"]), 
                    contig=contig, start=int(start), end=int(stop)
                )
                count_aln+=1
                if queryIDs:
                    if  attributes["Target"] in geneIDs_map_counts:
                        geneIDs_map_counts[attributes["Target"]].append(mini_alignment)
                    else:
                        geneIDs_map_counts[attributes["Target"]] = [mini_alignment]
                else:
                    geneIDs_map_counts[attributes["ID"]] = mini_alignment
        if queryIDs:
            print(f"read {count_aln} MINIPROT alignments from {len(geneIDs_map_counts)} query proteins")
        else:
            print(f"\t(read {count_aln} alignments from {miniprot_file})")

    return geneIDs_map_counts

if __name__ == "__main__":

    username = "miltr339"

    sex_chr_contigs = sex_chromosomes.get_contig_names()
    native_annot = get_gff_annot_paths(username=username)
    miniprot_annot = get_miniprot_annot_paths(username=username)
    single_exon_paths_dict_native,single_exon_paths_dict_miniprot = single_exon_paths(username=username)


    if True:
        ## plot single exon proportions:
        single_exon_dict_native, num_single_exon_dict_native, num_transcripts_dict_native = get_single_exon_genes(single_exon_paths_dict_native, write_to_file=False, outfile_name="native_single_exon_transcripts_list_14_species.txt", include_total_gene_num = True)
        single_exon_dict_miniprot, num_single_exon_dict_miniprot, num_transcripts_dict_miniprot = get_single_exon_genes(single_exon_paths_dict_miniprot, write_to_file=False, outfile_name="native_single_exon_transcripts_list_14_species.txt", include_total_gene_num = True)

        SE_numbers_native = {
            "num_single_exon_genes" : num_single_exon_dict_native,
            "num_total_genes" : num_transcripts_dict_native
        }
        SE_numbers_miniprot = {
            "num_single_exon_genes" : num_single_exon_dict_miniprot,
            "num_total_genes" : num_transcripts_dict_miniprot
        }
        print(f"\t --> native_numbers = {SE_numbers_native}\n\t --> miniprot_numbers = {SE_numbers_miniprot}\n")
        print()

        data_dir = f"/Users/{username}/work/PhD_code/PhD_chapter2/data/single_exon_stats"
        plot_single_exon_no_species_specific_three_annot(native_numbers=SE_numbers_native, miniprot_numbers=SE_numbers_miniprot, species_names=native_annot.keys(), filename = f"{data_dir}/single_exon_proportions.png")

    if False:
        counts_dict = {}
        for species, miniprot_file in miniprot_annot.items():
            print(f"=================== {species} ===================")
            counts_dict_species = miniprot_parse_alignment(miniprot_file = miniprot_file)
            keys = ["Aobt_anno1.g10014.t1"]# list(counts_dict_species.keys())
            print(f"{keys[0]} : {len(counts_dict_species[keys[0]])} alignments:")
            for aln in counts_dict_species[keys[0]]:
                print(aln)
            break 