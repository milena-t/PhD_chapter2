"""
This is code from the chapter III revision where I tested multi-species vs. pairwise dNdS. It makes 1-to-1 orthogroups from pairwise orthologs
I will find all other orhtologs of any A. obtectus gene, and then check if they are also each other's partners to identify groups of four 1-to-1 orthologs
Originally this requires exclusive X-linkage but I require at least one only
"""

import pandas as pd 

def get_paths(username = "miltr339"):
    dirname=f"/Users/{username}/work/chapter2/BRHs/"
    files_dict = {
        ("A_obtectus","B_siliquastri") : f"{dirname}A_obtectus_vs_B_siliquastri_BRHs.tsv",
        ("A_obtectus","B_varius") : f"{dirname}A_obtectus_vs_B_varius_BRHs.tsv",
        ("A_obtectus","C_chinensis") : f"{dirname}A_obtectus_vs_C_chinensis_BRHs.tsv",
        ("A_obtectus","C_maculatus") : f"{dirname}A_obtectus_vs_C_maculatus_BRHs.tsv",
        ("A_obtectus","D_carinulata") : f"{dirname}A_obtectus_vs_D_carinulata_BRHs.tsv",
        ("A_obtectus","D_sublineata") : f"{dirname}A_obtectus_vs_D_sublineata_BRHs.tsv",
        ("B_siliquastri","B_varius") : f"{dirname}B_siliquastri_vs_B_varius_BRHs.tsv",
        ("B_siliquastri","C_chinensis") : f"{dirname}B_siliquastri_vs_C_chinensis_BRHs.tsv",
        ("B_siliquastri","C_maculatus") : f"{dirname}B_siliquastri_vs_C_maculatus_BRHs.tsv",
        ("B_siliquastri","D_carinulata") : f"{dirname}B_siliquastri_vs_D_carinulata_BRHs.tsv",
        ("B_siliquastri","D_sublineata") : f"{dirname}B_siliquastri_vs_D_sublineata_BRHs.tsv",
        ("B_varius","C_chinensis") : f"{dirname}B_varius_vs_C_chinensis_BRHs.tsv",
        ("B_varius","C_maculatus") : f"{dirname}B_varius_vs_C_maculatus_BRHs.tsv",
        ("B_varius","D_carinulata") : f"{dirname}B_varius_vs_D_carinulata_BRHs.tsv",
        ("B_varius","D_sublineata") : f"{dirname}B_varius_vs_D_sublineata_BRHs.tsv",
        ("C_chinensis","C_maculatus") : f"{dirname}C_chinensis_vs_C_maculatus_BRHs.tsv",
        ("C_chinensis","D_carinulata") : f"{dirname}C_chinensis_vs_D_carinulata_BRHs.tsv",
        ("C_chinensis","D_sublineata") : f"{dirname}C_chinensis_vs_D_sublineata_BRHs.tsv",
        ("C_maculatus","D_carinulata") : f"{dirname}C_maculatus_vs_D_carinulata_BRHs.tsv",
        ("C_maculatus","D_sublineata") : f"{dirname}C_maculatus_vs_D_sublineata_BRHs.tsv",
        ("D_carinulata","D_sublineata") : f"{dirname}D_carinulata_vs_D_sublineata_BRHs.tsv",
    }
    return files_dict


class SpeciesOrtholog:
    """
    A species ortholog with a key of one geneID, species ident and a dict of orthologs { geneID : species }
    """
    def __init__(self, geneID:str, species:str, orthologs:dict, chromosome:str) -> None:
        self.geneID = geneID
        self.species=species
        self.orthologs=orthologs
        self.chromosome=chromosome

    def add_ortholog(self, key_value_tuple:tuple):
        orthologID,species=key_value_tuple
        self.orthologs[orthologID]= species

    def test_ortholog_presence(self, species:str, orthologID:str):
        if orthologID in self.orthologs:
            if self.orthologs[orthologID] == species:
                return True
            else:
                return False
        else:
            return False
    
    def __str__(self):
        return(
            f"""
gene {self.geneID} in species: {self.species} (on chromosome {self.chromosome})
\thas {len(self.orthologs)} orthologs: {self.orthologs}
            """
        )

class OrthoGroup:
    """
    group of orthologs that is 1-to-1 in all species as a dictionary like { species : geneID }
    as well as a dict of chromosome linkage { species : chr }
    """
    def __init__(self, groupID:int, species_dict:dict, chromosome:dict) -> None:
        self.groupID=groupID
        self.species_dict=species_dict
        self.chromosome=chromosome

    def __str__(self) -> str:
        return(f"""
1-to-1 orthogroup: {self.groupID}
is present in species {self.species_dict}
on chromosomes {self.chromosome}
        """)
    

def get_pair_partner(pair_tuple, species):
    pair_list = list(pair_tuple)
    pair_list.remove(species)
    assert len(pair_list)==1
    return pair_list[0]

def get_all_orthologs(pair_paths_dir, species_list, outfile, require_same_chr = False):
    """
    write a list of 1-to-1 orthologs across all species in species_list into an outfile
    """

    ## read all orthologs into class
    orthologs_dict = { species : {} for species in species_list}
    for species in species_list:
        species_orthologs = {}
        print(f"================= {species} =================")
        for pair_tuple, pair_filepath in pair_paths_dir.items():
            if species not in pair_tuple:
                continue

            species_partner = get_pair_partner(pair_tuple=pair_tuple, species=species)
            df=pd.read_csv(pair_filepath, sep="\t")

            chr_colname = "chromosome"
            if pair_tuple[1]==species:
                chr_colname = "chromosome.1"
            
            print(f"\t* {species_partner} : \t{pair_filepath}")
            for i,ortholog_line in df.iterrows():
                species_geneID = ortholog_line[species]
                partner_geneID = ortholog_line[species_partner]
                chromosome = ortholog_line[chr_colname]
                if species_geneID not in species_orthologs:
                    species_orthologs[species_geneID] = SpeciesOrtholog(geneID=species_geneID, species=species, orthologs={partner_geneID:species_partner}, chromosome=chromosome)
                else:
                    species_orthologs[species_geneID].add_ortholog(key_value_tuple = (partner_geneID,species_partner))

        orthologs_dict[species]=species_orthologs
        
    print("\n")
    for species, orthologs_dict_sp in orthologs_dict.items():
        print(f"{species} has {len(orthologs_dict_sp)} orthologs (Type: {type(orthologs_dict_sp)})")
    
    # loop through all orthologs in species_list[0] and check their orthologs with the full ortholog dictionary. 
    # Make list of all of these gene IDs and check if list(set()) gives the same length as species list (only four unique ones)
    # immediately skip past geneIDs that aren't represented in one other species
    # also skip past ones where the chromosome assignment is not uniquely A or X but mixed

    focal_species=species_list[0]
    num_species = len(species_list)
    skipped_focal_geneID_list = []
    wrong_chr = 0
    all_OGs = len(orthologs_dict[focal_species])

    one_to_ones = {}
    OG_id = 1
    for focal_geneID, ortholog_instance in orthologs_dict[focal_species].items():
        focal_chr = ortholog_instance.chromosome
        orthologs_list =[(focal_geneID,focal_species)] # make list of tuples (geneID,species) to avoid duplication overwrite stuff that will happen with dicts 

        # skip focal orthologs that don't have the correct number of partners
        if len(ortholog_instance.orthologs) != len(species_list)-1:
            skipped_focal_geneID_list.append(focal_geneID)
            continue
        
        chr_linkage_dict = {focal_species : focal_chr}
        
        # check the orthologs for all partners
        for partnerID,partner_species in ortholog_instance.orthologs.items():
            # this should not be possible, at least one reciprocal ortholog of focal_geneID should be in the partner species dict under partnerID
            if partnerID not in orthologs_dict[partner_species]:
                print(orthologs_dict[partner_species], len(orthologs_dict[partner_species]))
                print(f"focal species : {focal_species}")
                print(orthologs_dict[focal_species][focal_geneID])
                print(f"partner species : {partner_species}")
                print(orthologs_dict[partner_species][partnerID])
                raise RuntimeError(f"""
{focal_geneID} in {focal_species} does not have a reciprocal ortholog in {partner_species} for partner ortholog {partnerID}!
Something is wrong when making the dictionaries
                """)
            
            # skip partner orthologs that don't have the correct number of partners
            partner_orthologs_dict = orthologs_dict[partner_species][partnerID].orthologs
            if len(partner_orthologs_dict) != num_species-1:
                skipped_focal_geneID_list.append(focal_geneID)
                break  
            
            partner_chr = orthologs_dict[partner_species][partnerID].chromosome
            if require_same_chr:
                # skip partner orthologs that are not on the same chromosome as the focal ortholog
                if partner_chr != focal_chr:
                    skipped_focal_geneID_list.append(focal_geneID)
                    wrong_chr+=1
                    break
            else:
                chr_linkage_dict[partner_species] = partner_chr
            
            for geneID,species in partner_orthologs_dict.items():
                partner_tupel = (geneID,species)
                orthologs_list.append(partner_tupel)
        
        orthologs_list = list(set(orthologs_list))
        orthologs_species_list = list(set([species for gID,species in orthologs_list]))
        # skip when unique set of orthologs and species is not correct
        if len(orthologs_list) != num_species or len(orthologs_species_list) != num_species:
            # print(orthologs_list)
            skipped_focal_geneID_list.append(focal_geneID)
            continue
        
        # write to dict if all requirements are fulfilled
        onetoone_orthologs_dict = {species : geneID for geneID,species in orthologs_list}
        # assert len(chr_linkage_dict) == len(species_list)
        chr_linkage_dict = {species : orthologs_dict[species][geneID].chromosome for geneID,species in orthologs_list}
        one_to_ones[OG_id] = OrthoGroup(groupID=OG_id,species_dict=onetoone_orthologs_dict, chromosome=chr_linkage_dict)
        OG_id+=1

    skipped_focal_geneID_list = list(set(skipped_focal_geneID_list))
    intersection = list(set(skipped_focal_geneID_list) & set(one_to_ones))
    skipped_list_len = len(skipped_focal_geneID_list)
    sum_OGs = len(one_to_ones) + skipped_list_len
    if require_same_chr:
        wrongchr_str = f"({wrong_chr} were skipped because they were not exclusive to X or A in all species) "
    else:
        wrongchr_str = f"(no exclusive X or A linkage requirement)"
    print(f"""
{skipped_list_len} out of {all_OGs} A. obtectus orthologs do not fulfill the requirements for 1-to-1 in all species.
{wrongchr_str}
{len(one_to_ones)} 1-to-1 orthologs with all {len(species_list)} species detected
{len(one_to_ones)} + {skipped_list_len } = {sum_OGs} (should be {all_OGs}, difference is {sum_OGs - all_OGs}), length of intersection: {len(intersection)})
""")

    with open(outfile, "w") as output_file:
        header = ",".join(species_list)
        if require_same_chr:
            output_file.write(f"orthogroup_ID,chromosome,{header}\n") 
            for OG_id, orthogroup in one_to_ones.items():
                geneIDs = ",".join([orthogroup.species_dict[species] for species in species_list])
                outline = f"{orthogroup.groupID},{orthogroup.chromosome},{geneIDs}\n"
                output_file.write(outline)
        else:
            chr_link_header = ",".join([f"{species}:chr" for species in species_list])
            output_file.write(f"orthogroup_ID,{header},{chr_link_header}\n") 
            for OG_id, orthogroup in one_to_ones.items():
                geneIDs = ",".join([orthogroup.species_dict[species] for species in species_list])
                try:
                    chr_linkage = ",".join([orthogroup.chromosome[species] for species in species_list])
                except:
                    print(orthogroup)
                    raise RuntimeError(f"one species missing chr assignment")
                outline = f"{orthogroup.groupID},{geneIDs},{chr_linkage}\n"
                output_file.write(outline)



    print(f"\n-- written to output file: {outfile}")



if __name__ == "__main__":

    username="miltr339"
    brhs_tables = get_paths(username=username)

    get_all_orthologs(pair_paths_dir=brhs_tables, 
        species_list = ["A_obtectus","B_siliquastri","B_varius","C_chinensis","C_maculatus","D_carinulata","D_sublineata"], 
        outfile=f"/Users/{username}/work/PhD_code/PhD_chapter2/data/BRH_orthologs/chrysomelidae_orthologs.csv")