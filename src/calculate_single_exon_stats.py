### calculate gff statistics focused on single- and multi-exon genes
## write the results to stdout
"""
This is just a wrapper script around print_single_exon_stats() in parse_gff.py to make automating it a little easier
useage: 
python3 calculate_single_exon_stats.py path/to/annotation.gff [True|False]
True includes a list of the single-exon gene IDs in the output, False gives just the statistics
"""

import parse_gff as gff
import sys

if len(sys.argv)>1:
    print(f"{sys.argv}\n\n")
    if len(sys.argv)==3:
        annot_filepath = sys.argv[1]
        arg =  sys.argv[2]
        topfeat = "false"
    elif len(sys.argv)==4:
        annot_filepath = sys.argv[1]
        arg =  sys.argv[2]
        topfeat = sys.argv[3]
    else:
        print(f"\narguments: {sys.argv}")
        print(f"useage: python3 calculate_single_exon_stats.py annotation.gff [True|False]")
        sys.exit(1)

    if arg.lower() == "true":
        boolean_value = True
    elif arg.lower() == "false":
        boolean_value = False
    else:
        print("Invalid second argument. Please pass 'True' or 'False'.")
        sys.exit(1)

    if topfeat.lower() == "true":
        boolean_topfeat = True
        print(f"*** top-level feature mRNA")
    elif topfeat.lower() == "false":
        boolean_topfeat = False
    else:
        print("Invalid third argument. Please pass 'True' or 'False'.")
        sys.exit(1)


    filename = annot_filepath.split("/")[-1]
    print(f"{filename}")

    gff.print_single_exon_stats(filepath=annot_filepath, include_list=boolean_value, mRNA_top_feature=boolean_topfeat) # include or exclude a list of all the single-exon IDs
else:
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
    
    for species,annot_filepath in get_miniprot_annot_paths().items():
        print(f"======================= {species} =======================")
        gff.print_single_exon_stats(filepath=annot_filepath, include_list=True, mRNA_top_feature=True) # include or exclude a list of all the single-exon IDs
        