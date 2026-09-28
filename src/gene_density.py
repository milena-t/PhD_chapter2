"""
Calculate gene density (via exon annotated regions)
"""

import sex_chromosomes

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
        "A_obtectus" : f"{filedir}A_obtectus.masked_align.gff",
        "B_siliquastri" : f"{filedir}B_siliquastri.masked_align.gff",
        "B_varius" : f"{filedir}B_varius.masked_align.gff",
        "C_chinensis" : f"{filedir}C_chinensis.masked_align.gff",
        "C_maculatus" : f"{filedir}C_maculatus.masked_align.gff",
        "D_carinulata" : f"{filedir}D_carinulata.masked_align.gff",
        "D_sublineata" : f"{filedir}D_sublineata.masked_align.gff",
    }
    return outdict


if __name__ == "__main__":

    username = "miltr339"

    sex_chr_contigs = sex_chromosomes.get_contig_names()
    native_annot = get_gff_annot_paths(username=username)
