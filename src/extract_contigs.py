"""
Extract all contigs from a multifasta assembly into a separate fasta file.

Extract X and Y contigs from Cmac Y_s and Y_l assemblies to make separate alignments
"""

from Bio import SeqIO
import sex_chromosomes


sex_chr_contigs = sex_chromosomes.Cmac_S_L_nonscaffolded_contig_names()

ass_dir = "/proj/coleoptera-genomics-2025/snic2021-6-30/Philipp/HiFi/Whole_genome_alinment/Alinment_pt_036"
assemblies = {
    "Y_s" : f"{ass_dir}/pt_036_001_hifiasm_primary.fasta",
    "Y_l" : f"{ass_dir}/pt_036_002_hifiasm_primary.fasta",
}

def extract_contigs(main_fasta:str, contig_list:list, out_fasta:str):
    outfile_records = []
    for record in SeqIO.parse(main_fasta, "fasta"):
        if record.id in contig_list:
            outfile_records.append(record)
    SeqIO.write(outfile_records, out_fasta, "fasta")



if __name__ == "__main__":

    for haplotype in ["Y_s","Y_l"]:
        assembly = assemblies[haplotype]