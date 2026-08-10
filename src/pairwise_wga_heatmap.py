"""
Plot a heatmap for all pairwise genome alignments. I will be inspired by this: https://www.nature.com/articles/s41586-024-07473-2/figures/1
and they do: "percentage of reference, as shown on the x axis, covered by the query, as shown on the y axis"
"""

import sex_chromosomes
import pandas as pd

def fasta_indices(username="milena"):
    """
    read the assembly index files and make dicts of {
        species : {contig : length, contig : length , ...}
    }
    rsync -azP "milenatr@pelle.uppmax.uu.se:/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/assemblies/*fai" /Users/miltr339/work/chapter2/wga/
    """
    dirname = f"/Users/{username}/work/chapter2/wga/"
    filenames_dict = {
        "A_obtectus" : f"{dirname}A_obtectus.masked.fna.fai",
        "B_siliquastri" : f"{dirname}B_siliquastri.masked.fna.fai",
        "B_varius" : f"{dirname}B_varius.masked.fna.fai",
        "C_chinensis" : f"{dirname}C_chinensis.masked.fna.fai",
        "C_maculatus" : f"{dirname}C_maculatus.masked.fna.fai",
        "D_carinulata" : f"{dirname}D_carinulata.masked.fna.fai",
        "D_sublineata" : f"{dirname}D_sublineata.masked.fna.fai",
    }
    def read_faidx(infile_path):
        df = pd.read_csv(infile_path, sep="\t", header=None)
        return {contig : length for contig,length in zip(df.iloc[:,0],df.iloc[:,1])}

    out_dict = {species : read_faidx(filename) for species,filename in filenames_dict.items()}
    return out_dict

def get_aln_coord_files(username="milena"):
    dirname = f"/Users/{username}/work/chapter2/wga/"
    out_dict = {
        "A_obtectus" : {
            "B_siliquastri" : f"{dirname}ref_A_obtectus.masked_query_B_siliquastri.masked.delta_filtered_coords",
            "B_varius" : f"{dirname}",
            "C_chinensis" : f"{dirname}",
            "C_maculatus" : f"{dirname}",
            "D_carinulata" : f"{dirname}",
            "D_sublineata" : f"{dirname}",
        },
        "B_siliquastri" : {
            "A_obtectus" : f"{dirname}",
            "B_varius" : f"{dirname}",
            "C_chinensis" : f"{dirname}",
            "C_maculatus" : f"{dirname}",
            "D_carinulata" : f"{dirname}",
            "D_sublineata" : f"{dirname}",
        },
        "B_varius" : {
            "A_obtectus" : f"{dirname}",
            "B_siliquastri" : f"{dirname}",
            "C_chinensis" : f"{dirname}",
            "C_maculatus" : f"{dirname}",
            "D_carinulata" : f"{dirname}",
            "D_sublineata" : f"{dirname}",

        },
        "C_chinensis" : {
            "A_obtectus" : f"{dirname}",
            "B_siliquastri" : f"{dirname}",
            "B_varius" : f"{dirname}",
            "C_maculatus" : f"{dirname}",
            "D_carinulata" : f"{dirname}",
            "D_sublineata" : f"{dirname}",
        },
        "C_maculatus" : {
            "A_obtectus" : f"{dirname}",
            "B_siliquastri" : f"{dirname}",
            "B_varius" : f"{dirname}",
            "C_chinensis" : f"{dirname}",
            "D_carinulata" : f"{dirname}",
            "D_sublineata" : f"{dirname}",
        },
        "D_carinulata" : {
            "A_obtectus" : f"{dirname}",
            "B_siliquastri" : f"{dirname}",
            "B_varius" : f"{dirname}",
            "C_chinensis" : f"{dirname}",
            "C_maculatus" : f"{dirname}",
            "D_sublineata" : f"{dirname}",
        },
        "D_sublineata" : {
            "A_obtectus" : f"{dirname}",
            "B_siliquastri" : f"{dirname}",
            "B_varius" : f"{dirname}",
            "C_chinensis" : f"{dirname}",
            "C_maculatus" : f"{dirname}",
            "D_carinulata" : f"{dirname}",
        }
    }
    return out_dict

def count_ref_aln_coordinates(aln_filtered_coords_path, ref_faidx_dict):
    """
    get the percentage of bases of each reference contig that are covered by the query. 
    full contig length is determined from the faidx file
    The first two column in the alignment coordinates are start and end of the alignment in the reference contig
    """
    out_cov_bases = {}
    with open(aln_filtered_coords_path, "r") as aln_filtered_coords:
        for line in aln_filtered_coords.readlines():
            refstart,refend,qstart,qend,reflen,qlen,percident,refcontig,qcontig = line.strip().split("\t")                
            if refcontig in out_cov_bases:
                out_cov_bases[refcontig] += int(reflen)
            else:
                out_cov_bases[refcontig] = int(reflen)

    out_perc_bases = {contig : num_bases / ref_faidx_dict[contig] for contig,num_bases in out_cov_bases.items()}
    return out_perc_bases
    

if __name__ == "__main__":
    username = "milena"
    sex_chr_dict = sex_chromosomes.get_contig_names()
    faidx_dicts = fasta_indices(username=username)
    aln_coord_files = get_aln_coord_files(username=username)
    test_perc_bases = count_ref_aln_coordinates(aln_filtered_coords_path=aln_coord_files["A_obtectus"]["B_siliquastri"], ref_faidx_dict=faidx_dicts["A_obtectus"])
    for contig,perc in test_perc_bases.items():
        print(f"{contig} : {perc:.4f}")

