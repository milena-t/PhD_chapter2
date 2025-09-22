"""
MCScanX doesn't find any collinear genes and I suspect it has something to do with the contig names in my assemblies
The github documentation specifically says they should be e.g. for A. thaliana chromosome 1:  'at1' 

Therefore I will modify the contig names accordingly for all species and generate a lookup table to backtrace them
"""


def make_lookup_table(annotation_path:str, species_initials:str, lookup_table_name:str = ""):
    """
    generate a lookup table for contig names from an annotation to species_index and a number
    """
    if lookup_table_name == "":
        annotation_name = annotation_path.split("/")[-1]
        lookup_table_name = f"{annotation_name}.MCScanX_contig_name_lookup_table.tsv"

    ## get list of all contigs
    with open(annotation_path, "r") as annotation_file:
        # read file and remove comment lines
        annotation_lines = [line for line in annotation_file.readlines() if line[0]!="#"]
        contig_names = annotation_lines
        for i,annotation_line in enumerate(annotation_lines):
            line=annotation_line.strip().split("\t")
            contig_names[i] = line[0]
        
        contig_names = list(set(contig_names))
    
    # write to lookup table file
    with open(lookup_table_name, "w") as lookup_table_file:
        for i, contig_name in enumerate(contig_names):
            mcscanx_name = f"{species_initials}{str(i+1 )}"
            lookup_line = f"{contig_name}\t{mcscanx_name}\n"
            lookup_table_file.write(lookup_line)
    
    print(f"lookup table with species initial {species_initials} written to {lookup_table_name}")
    return lookup_table_name


def read_lookup_table_to_dict(lookup_dict_path:str):
    """
    read lookup table to dict { original_contig : new_contig }
    """
    lookup_dict = {}
    with open(lookup_dict_path, "r") as lookup_dict_file:
        for line in lookup_dict_file.readlines():
            original_contig, new_contig = line.strip().split("\t")
            lookup_dict[original_contig] = new_contig
    return lookup_dict


def modify_bedfile_with_lookup_table(lookup_dict_path:str, bedfile_path:str):
    """ 
    modify the bedfile input to contain the updated contig names
    """
    bedfile_new_path = bedfile_path[:-4]
    bedfile_new_path = f"{bedfile_new_path}_uniform_contig_names.bed"
    lookup_dict = read_lookup_table_to_dict(lookup_dict_path=lookup_dict_path)

    with open(bedfile_path, "r") as bedfile_old, open(bedfile_new_path, "w") as bedfile_new:
        for bed_line in bedfile_old.readlines():
            contig, start, stop, gene = bed_line.strip().split("\t")
            contig_new = lookup_dict[contig]
            bedfile_new_line = f"{contig_new}\t{gene}\t{start}\t{stop}\n"
            bedfile_new.write(bedfile_new_line)
    print(f"updated bedfile written to {bedfile_new_path}")



if __name__ == "__main__":
    
    annotation_file = "/Users/miltr339/work/native_annotations/a_obtectus_native_isoform_filtered.gff"
    bedfile = "/Users/miltr339/work/native_annotations/a_obtectus_native_isoform_filtered_simplified.bed"
    lookup_table = make_lookup_table(annotation_path=annotation_file, species_initials="ao")

    modify_bedfile_with_lookup_table(lookup_dict_path=lookup_table, bedfile_path=bedfile)

