import numpy as np

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