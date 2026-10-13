"""
make fasta files with the two BRH sequences from the blast results
to calculate the dNdS ratio
"""

from Bio import SeqIO
import pandas as pd
from tqdm import tqdm
import os 

def nucleotides_paths(username = "miltr339"):
    nuc_dir = f"/Users/{username}/work/chapter2/nucleotide_data/"
    filepaths_dict = {
        "A_obtectus" : f"{nuc_dir}A_obtectus.fna",
        "B_varius" : f"{nuc_dir}B_varius.fna",
        "B_siliquastri" : f"{nuc_dir}B_siliquastri.fna",
        "C_chinensis" : f"{nuc_dir}C_chinensis.fna",
        "C_maculatus" : f"{nuc_dir}C_maculatus.fna",
        "D_carinulata" : f"{nuc_dir}D_carinulata.fna",
        "D_sublineata" : f"{nuc_dir}D_sublineata.fna",
    }
    return filepaths_dict



def make_bruchini_ortholog_fasta_files(brh_path, nucleotides_dict, species_list, outdir = "", chr_type="A"):
    """
    make fasta files for dNdS from pairwise blast results
    with chromosome type you can restrict the orthologs to X-exclusive or A-exclusive by setting "X" or "A"
    """

    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    ### make dictionary of { species : { transcriptID : SeqIO_record } }
    ### to quickly access all the sequences when making fasta files
    print(f"\nmaking nucleotide fasta files for orthologs")
    nucleotides_seqs_dict = {}
    print(f"reading nucleotide files...")
    for species in species_list:
        # make dict from the nucleotide files
        if species not in nucleotides_seqs_dict.keys():
            nucleotides_seqs_dict[species] = {}
            for cds_record in SeqIO.parse(nucleotides_dict[species], "fasta"):
                nucleotides_seqs_dict[species][cds_record.id] = cds_record


    brh_table = pd.read_csv(brh_path, sep=",")
    
    pair_dirname = f"{outdir}"
    if not os.path.isdir(pair_dirname):
        os.mkdir(pair_dirname)
    
    not_found_partners_dict = {species : 0 for species in species_list}
    num_species = len(species_list)
    count_pairs = 0
    for i,brh_set in tqdm(brh_table.iterrows()):
        # print(f"ID1 = {brh_pair[1]}, ID2 = {brh_pair[3]}")
        brh_seq_records = []

        for species in species_list:
            try:
                record1 = nucleotides_seqs_dict[species][brh_set[species]]
                record1.id=f"{brh_set[species]}_{species}"
                brh_seq_records.append(record1)
            except:
                not_found_partners_dict[species]+=1
                print(brh_set)
                try:
                    transcripts_fasta = nucleotides_dict[species]
                except:
                    raise RuntimeError(f"{species} does not have a transcripts file in 'nucleotides_dict'.")
                # raise RuntimeError(f"{brh_pair[1]} not found in {transcripts_fasta}")
            
            id = brh_set["orthogroup_ID"]
        if len(brh_seq_records)==num_species:
            fasta_name = f"{pair_dirname}Chrysomelidae_bruchini_{chr_type}-linked_ortholog_{id}.fasta"
            SeqIO.write(sequences=brh_seq_records, handle=fasta_name, format="fasta")
    
            count_pairs += 1

    # print(f"\t* {species_partner}: {count_pairs} orthologs")
    print(f"\t\t{count_pairs} counted, {not_found_partners_dict} transcripts not found")
    print(f"\ndone! all outfiles written to \n{outdir}")






if __name__ == "__main__":
    username = "miltr339"
    nucleotides_dict = nucleotides_paths(username)
    
    if True:
        diorhabda = ["D_carinulata","D_sublineata"]
        bruchini = ["C_maculatus","C_chinensis","A_obtectus","B_siliquastri","B_varius"]

        bruchini_A = f"/Users/{username}/work/PhD_code/PhD_chapter2/data/BRH_orthologs/chrysomelidae_orthologs_A_in_bruchini.csv"
        bruchini_X = f"/Users/{username}/work/PhD_code/PhD_chapter2/data/BRH_orthologs/chrysomelidae_orthologs_X_in_bruchini.csv"
        species_list = diorhabda+bruchini
        outdir_X = f"/Users/{username}/work/chapter2/BRHs/multifasta_diorhabdaX_bruchiniX/"
        outdir_A = f"/Users/{username}/work/chapter2/BRHs/multifasta_diorhabdaX_bruchiniA/"
        make_bruchini_ortholog_fasta_files(brh_path=bruchini_A, nucleotides_dict=nucleotides_dict, species_list=species_list, chr_type="A", outdir= outdir_A)
        make_bruchini_ortholog_fasta_files(brh_path=bruchini_X, nucleotides_dict=nucleotides_dict, species_list=species_list, chr_type="X", outdir= outdir_X)
    
    