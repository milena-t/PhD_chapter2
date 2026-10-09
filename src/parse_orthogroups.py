import numpy as np
import pandas as pd
from statistics import mean
import linkage_groups as lg
import parse_gff as gff
import miniprot_stats_comparison as minialn
from collections import defaultdict


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


def parse_orthogroups_dict(filepath, species_list ="", OG_header = "HOG"):
    orthogroups_df = pd.read_csv(filepath, sep="\t")
    headers = set(orthogroups_df.columns)
    if species_list != "":
        annot_species = set(species_list)
    else:
        annot_species = set(list(headers)[3:])
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
    def __init__(self, transcript_ID:str, species:str, LG:lg.LinkageGroup, is_miniprot=False, gff_feature:gff.Feature=None) -> None:
        self.transcript_ID=transcript_ID
        self.species=species
        self.LG=LG
        self.is_miniprot=is_miniprot
        self.gff_feature=gff_feature
    
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

    def miniprot_members(self):
        return [m for m in self.members if m.is_miniprot]

    def sex_linkage(self, species = ""):
        """
        return a dict of geneID lists that are X/Y or O linked (O for other, either A or unplaced)
        """
        if species !="":
            out_dict = {lg_ : 0 for lg_ in self.linkagegroups(species)}
            for m in self.members:
                if m.species==species:
                    out_dict[m.LG]+=1
            return out_dict
        else:
            out_dict ={ s: {"X" : 0, "Y" : 0, "O" : 0} for s in self.species()}
            for m in self.members:
                if m.LG == lg.LinkageGroup.LGX:
                    out_dict[m.species]["X"] += 1
                if m.LG == lg.LinkageGroup.LGY:
                    out_dict[m.species]["Y"] += 1
                if m.LG != lg.LinkageGroup.LGY and m.LG != lg.LinkageGroup.LGX:
                    out_dict[m.species]["O"] += 1
            return out_dict


    def linkagegroups(self, species=""):
        if species =="":
            return list(set([m.LG for m in self.members]))
        else:
            return list(set([m.LG for m in self.members if m.species==species]))

    def species(self):
        return list(set([m.species for m in self.members]))

    def __str__(self) -> str:
        outlink = {s: [] for s in self.species()}
        for m in self.members:
            outlink[m.species].append(m.LG)
        outstr = "\n  * ".join([f"{s} ({len(lglist)}) : {set(lglist)}" for s,lglist in outlink.items()])

        minilink = {s: [] for s in self.species()}
        for m in self.miniprot_members():
            minilink[m.species].append(m.gff_feature.ID)
        ministr = "\n  - ".join([f"{s} ({len(lglist)}) : {set(lglist)}" for s,lglist in minilink.items()])
        return (
f"""
Orthogroup {self.OG_id} has {self.OG_size()} members and is present on {len(self.species())} species ({self.species()})
it is on {len(self.linkagegroups())} linkage groups ({self.linkagegroups()})
in these species:
  * {outstr}
with {len(self.miniprot_members())} being miniprot alignments:\n  - {ministr}
"""
        )
    def __repr__(self):
        return "OrthoGroup"

    

def parse_orthogroups_class(filepath, annot_species, annotations_dict, unassigned_genes_path = "", miniprot_paths_dict = {}, add_gff_feature=False):
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
        miniprot_dict = { species : minialn.miniprot_parse_alignment(miniprot_paths_dict[species], nest_cds=True) for species in annot_species}
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
                if add_gff_feature:
                    geneIDs_species = [ 
                        OGMember(
                            transcript_ID=gid, 
                            species=species, 
                            LG=lg.assign_linkagegoup(annot_gff_dict[species][gid].contig), 
                            gff_feature=annot_gff_dict[species][gid],
                            is_miniprot=False
                            ) 
                        for gid in geneIDs]
                else:
                    geneIDs_species = [ 
                        OGMember(
                            transcript_ID=gid, 
                            species=species, 
                            LG=lg.assign_linkagegoup(annot_gff_dict[species][gid].contig),
                            is_miniprot=False) 
                        for gid in geneIDs]
                
                all_geneids_class.extend(geneIDs_species)
        
        og_class_dict[HOG_id] = OrthoGroup(OG_id=HOG_id, members=all_geneids_class)

        if miniprot_paths_dict != {}:
            for member in og_class_dict[HOG_id].members:
                geneID = member.transcript_ID
                mini_species = member.species

                if geneID in miniprot_dict[mini_species]:
                    for miniprot_ID in miniprot_dict[mini_species][geneID]:
                        mini_contig = miniprot_ID.contig
                        og_class_dict[HOG_id].add_member(OGMember(
                            transcript_ID=miniprot_ID, 
                            species=mini_species, 
                            LG=lg.assign_linkagegoup(mini_contig), 
                            is_miniprot=True,
                            gff_feature=miniprot_ID,
                            ),
                        )
    
    return(og_class_dict)


class EarliestDuplication:
    def __init__(self, HOG_ID, geneIDs_side1:list, geneIDs_side2:list, support:float, tree_node:str=None) -> None:
        self.HOG_ID=HOG_ID
        self.geneIDs_side1=geneIDs_side1
        self.geneIDs_side2=geneIDs_side2
        self.support=support
        self.tree_node=tree_node



def parse_duplications(duplications_path, orthogroups_path, min_support = 0.5 ):
    
    ### make lookup table for {geneID : orthogroup}
    orthogroup_hogdict = parse_orthogroups_dict(orthogroups_path)
    # print(orthogroup_hogdict["N0.HOG0000001"])
    geneid_hogid_dict = {}
    for HOG_id, species_dict in orthogroup_hogdict.items():
        for species, geneIDs in species_dict.items():
            if isinstance(geneIDs, str): ## if not np.nan
                geneIDs_ = geneIDs.strip().split(", ")
                for geneID_ in geneIDs_:
                    sp_geneID = f"{species}_{geneID_}"
                    geneid_hogid_dict[sp_geneID] = HOG_id

    # print(geneid_hogid_dict["B_siliquastri_BRAKERILHT00000009893"])
    ## read duplications, make depth metric (lower is closer to root and therefore earlier) from node names
    df = pd.read_csv(duplications_path, sep="\t")
    df = df[df["Support"]>min_support]
    df = df[df["Type"] != "Terminal"] # remove terminal duplications since i don't care about those for sure
    df["Depth"] = df["Species Tree Node"].apply(lambda x : int(x[-1]) if len(x)==2 else x)
    df["Genes 1"] = df["Genes 1"].apply(lambda x : x.split(", "))
    df["Genes 2"] = df["Genes 2"].apply(lambda x : x.split(", "))

    ## associate HOG ideas with geneIDs
    def groug_by_HOG(geneIDs_list):
        d = defaultdict(list)
        for g in geneIDs_list:
            d[geneid_hogid_dict.get(g, "unassigned")].append(g)
        return dict(d)
    df["HOG 1"] = df["Genes 1"].apply(groug_by_HOG)
    df["HOG 2"] = df["Genes 2"].apply(groug_by_HOG)

    ## get unassigned genes
    unassigned_genes = set()
    for col in ("HOG 1", "HOG 2"):
        for d in df[col]:
            unassigned_genes.update(d.get("unassigned", []))
    unassigned_dict = defaultdict(list)
    for gid in unassigned_genes:
        species = gff.split_at_second_occurrence(gid)
        gid_ = gid.replace(f"{species}_", "")
        unassigned_dict.setdefault(species, [gid_]).append(gid_)
    unassigned_dict = {species : list(set(l)) for species,l in unassigned_dict.items()}

    print(f"---------------------------------")
    print(f"...reading {duplications_path}")
    for species,ulist in unassigned_dict.items():
        print(f"{species} ({len(ulist)} unassigned to a HOG)")
    print(f"---------------------------------")
    

    ## filter within HOG for the earliest HOG duplication eventd (oldest duplication on the gene tree, lowest node integer in gt_num)
    df["gt_num"] = df["Gene Tree Node"].str.extract(r"(\d+)", expand=False).astype(int) 
    # all unique HOGs involved in either side of the duplicaiton event
    df["HOG"] = df.apply(lambda r: set(r["HOG 1"]) | set(r["HOG 2"]), axis=1)

    long = df.explode("HOG")
    long = long[long["HOG"] != "unassigned"]
    earliest = (long.sort_values(["HOG", "Depth", "gt_num", "Support"],
                         ascending=[True, True, True, False]).drop_duplicates("HOG", keep="first").copy())

    ### make new dict with earliest duplication events for every HOG
    hog_events = defaultdict(list)
    hog_read=[]
    for r in earliest.to_dict("records"): 
        hog = r["HOG"]
        if hog not in hog_read:
            hog_events[hog] = EarliestDuplication(
                HOG_ID=hog, 
                geneIDs_side1=r["HOG 1"].get(hog, []),
                geneIDs_side2=r["HOG 2"].get(hog, []),
                support=r["Support"],
                tree_node=r["Species Tree Node"],
            )
        else:
            e = earliest[earliest["HOG"]==hog]
            raise RuntimeError(f"{hog} is duplicated even after filtering!\n{e}")
        
    return hog_events, unassigned_dict