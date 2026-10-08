import numpy as np
import pandas as pd
from statistics import mean
import linkage_groups as lg
import parse_gff as gff
import miniprot_stats_comparison as minialn

def get_orthogroup_sizes(orthogroup_dict, q = 0):
    """
    returns a dict of all orthogroups sizes. Do NOT confuse it with gene family sizes! (get_GF_sizes function instead)
    if a percentile q is specified other than 0, it returns only orthogroups whose size is > than the pth percentile of the size distribution
    """
    
    OG_sizes_dict = {}
    sizes = []
    for OG_id, species_dict in orthogroup_dict.items():
        size = 0
        for transcripts_list in species_dict.values():
            size += len(transcripts_list)
        
        OG_sizes_dict[OG_id] = size
        sizes.append(size)
    
    if q == 0:
        return OG_sizes_dict
    else:
        OG_sizes_filtered = {}
        sizes = np.array(sizes)
        percentile_size = np.percentile(sizes, q = q)
        for OG_id, size in OG_sizes_dict.items():
            if size > percentile_size:
                OG_sizes_filtered[OG_id] = size

        return OG_sizes_filtered


def parse_orthogroups_dict(filepath, species_list, OG_header = "HOG"):
    orthogroups_df = pd.read_csv(filepath, sep="\t")
    headers = set(orthogroups_df.columns)
    annot_species = set(species_list)
    headers_keep = list(annot_species & headers)
    headers = [f"{OG_header}"] + headers_keep
    
    df = orthogroups_df[headers]
    df = df.set_index(f"{OG_header}", drop=True)
    df_dict = df.to_dict(orient="index")
    # if OG_header != "HOG":
    #     print(df.loc["OG0014394"])
    
    return df_dict

def get_GF_sizes(orthogroups_dict):
    """
    modify the dictionary from the parse_orthogroups function
    
    {
        orthogroup : {
            species1 : num,
            species2 : num,
            ...
        }
        orthogroup2 : {
            species1 : num,
            species2 : num,
            ...
        }
    }
    """
    orthogroups_modified = {}
    for orthogroup_id , species_dict in orthogroups_dict.items():
        species_dict_nums = { species : len(transcripts_list) for species, transcripts_list in species_dict.items() if transcripts_list!=['']}
        orthogroups_modified[orthogroup_id] = species_dict_nums

    return orthogroups_modified

            
def get_mean_GF_size(orthogroups_dict:dict, species:str):
    GF_sizes = []
    for gene_families in orthogroups_dict.values():
        try:
            num_families = len(gene_families[species])
        except:
            num_families = 0
        GF_sizes.append(num_families)
    
    return mean(GF_sizes)
    
    
def get_species_in_OG_dict(OG_dict:dict) -> list:
    """
    get the list of all species included in the orthofinder dict
    """
    all_species = []
    for OG_id, GF_dict in OG_dict.items():
        species = list(GF_dict.keys())
        all_species.extend(species)
    all_species = list(set(all_species))
    return(all_species)


class OGMember:
    def __init__(self, transcript_ID:str, species:str, LG:lg.LinkageGroup, is_miniprot=False) -> None:
        self.transcript_ID=transcript_ID
        self.species=species
        self.LG=LG
        self.is_miniprot=is_miniprot
    
    def __repr__(self):
        return "OrthoGroupMember"



class OrthoGroup:
    def __init__(self, OG_id:str, members:list[OGMember]) -> None:
        self.OG_id=OG_id
        self.members=members

    def add_member(self, member:OGMember):
        self.members.append(member)

    def OG_size(self):
        return len(self.members)

    def GF_size(self,species=""):
        if species =="":
            sizes_dict = {s : 0 for s in self.species()}
            for m in self.members:
                sizes_dict[m.species] += 1
            return sizes_dict
        else:
            gfsize = 0
            for m in self.members:
                if m.species==species:
                    gfsize += 1
            return gfsize

    def linkagegroups(self):
        return list(set([m.LG for m in self.members]))

    def species(self):
        return list(set([m.species for m in self.members]))

    def __str__(self) -> str:
        outlink = {s: [] for s in self.species()}
        for m in self.members:
            outlink[m.species].append(m.LG)
        outstr = "\n  * ".join([f"{s} ({len(lglist)}) : {set(lglist)}" for s,lglist in outlink.items()])

        return (
f"""
Orthogroup {self.OG_id} has {self.OG_size()} members and is present on {len(self.species())} species ({self.species()})
it is on {len(self.linkagegroups())} linkage groups ({self.linkagegroups()})
in these species:
  * {outstr}\n
"""
        )
    def __repr__(self):
        return "OrthoGroup"

    

def parse_orthogroups_class(filepath, annot_species, annotations_dict, unassigned_genes_path = "", miniprot_paths_dict = {}):
    """
    parse orthogroups into the class structure with defined linkage group membership
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


    if miniprot_paths_dict != {}:
        print(f"--- read miniprot annotations ---")
        miniprot_dict = { species : minialn.miniprot_parse_alignment(miniprot_paths_dict[species]) for species in annot_species}
        print(f"---------------------------------")
    else:
        miniprot_dict = { species : [] for species in annot_species}
    
    og_class_dict = {}

    for HOG_id , OG_dict in annot_HOG_dict.items():
        sexchr_counts_dict = {}
        # print(f" --- {HOG_id} --- ")
        all_geneids_class = []
        for species, geneIDs in OG_dict.items():
            if isinstance(geneIDs, str): ## if not np.nan
                geneIDs = geneIDs.strip().split(", ")
                # print(f"   - {species} : {len(geneIDs)}")
                geneIDs_species = [ OGMember(transcript_ID=gid, species=species, LG=lg.assign_linkagegoup(annot_gff_dict[species][gid].contig)) for gid in geneIDs]
                all_geneids_class.extend(geneIDs_species)
        
        og_class_dict[HOG_id] = OrthoGroup(OG_id=HOG_id, members=all_geneids_class)

        if miniprot_paths_dict != {}:
            for member in og_class_dict[HOG_id].members:
                geneID = member.transcript_ID
                mini_species = member.species

                if geneID in miniprot_dict[mini_species]:
                    for miniprot_ID in miniprot_dict[mini_species][geneID]:
                        mini_contig = miniprot_ID.contig
                        og_class_dict[HOG_id].add_member(OGMember(transcript_ID=miniprot_ID, species=mini_species, LG=lg.assign_linkagegoup(mini_contig), is_miniprot=True))
    
    return(og_class_dict)