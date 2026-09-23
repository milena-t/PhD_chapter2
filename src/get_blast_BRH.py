"""
Get the transcript IDs of best reciprocal hits (BRH) of two species. the input is the outfmt 6 blast results in either direction
useage: python3 get_blast_BRH.py species1_species2.blast species2_species1.blast
and the output is a tsv file with two columns, one is species1 and one is species2, and each row contains a set of gene IDs that are each others BRH
I identify the best hit by the highest bit score
"""
from sex_chromosomes import get_contig_names as sex_chromosome_names
from circos.make_circos_karyotype_file import autosomes_lists
import parse_gff as gff
import argparse

def parse_args():
    # Create the parser
    program_description = """
script to compute best reciprocal hits and their chromosome location from a bidirectional all-vs-all proteinblast search and the corresponding species annotations
the result is a tsv file that includes the transcript IDs of the best hits and their chromosome location as X/Y/A (not actual contig IDs!)
"""
    parser = argparse.ArgumentParser(description=program_description)

    # Add the arguments
    parser.add_argument('--blast1', type=str, required=True, help='Absolute filepath to the all-vs-all blastp output where species1 is the query (outfmt 6')
    parser.add_argument('--blast2', type=str, required=True, help='Absolute filepath to the all-vs-all blastp output where species2 is the query (outfmt 6')
    parser.add_argument('-o', '--outfile', type=str, help='output filename, default: filename from blast1 + _BRH.tsv')
    parser.add_argument('--annotation1', type=str, required=True, help='Absolute filepath to the annotation that blast1 query is based on')
    parser.add_argument('--annotation2', type=str, required=True, help='Absolute filepath to the annotation that blast2 query is based on')
    parser.add_argument('--verbose', action='store_true', help="enable verbose mode")
    parser.add_argument('--X_contigs1', type=str, help='comma separated list of the names of X-linked contigs in species 1 (given directly in the command line with no spaces, not a path to a file!')
    parser.add_argument('--X_contigs2', type=str, help='comma separated list of the names of X-linked contigs in species 2 (given directly in the command line with no spaces, not a path to a file!')
    parser.add_argument('--Y_contigs1', type=str, help='comma separated list of the names of Y-linked contigs in species 1 (given directly in the command line with no spaces, not a path to a file!')
    parser.add_argument('--Y_contigs2', type=str, help='comma separated list of the names of Y-linked contigs in species 2 (given directly in the command line with no spaces, not a path to a file!')
 

    # Parse the arguments
    args = parser.parse_args()

    # set default values for non-obligatory arguments
    if not args.outfile:
        args.outfile = f"{args.blast1}_BRH.tsv"

    return args

blast_outfmt6_headers = ["qseqid", "rseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]


class BestHit:
    """
    which reference ID is the best hit to the query ID
    """
    
    def __init__(self, qseqid:str, rseqid:str, bitscore:float) -> None:
        self.qseqid = qseqid
        self.rseqid = rseqid
        self.bitscore = bitscore

    def update_besthit(self, new_qseqid:str, new_rseqid:str, new_bitscore:float) -> bool:
        assert new_qseqid == self.qseqid
        if new_bitscore > self.bitscore:
            self.rseqid = new_rseqid
            self.bitscore = new_bitscore
            return True
        else:
            return False

    def __str__(self) -> str:
        return(f"  * {self.qseqid} has the best hit {self.rseqid} with a bitscore of {self.bitscore}")


def read_best_hits(blast_infile_path:str) -> dict:
    """
    make a dictionary of all the best hits instances
    { queryID : BestHit ,  ...}
    """
    best_hits_dict = {}
    with open(blast_infile_path, "r") as blast_infile:
        blast_lines = blast_infile.readlines()
        for line in blast_lines:
            line = line.strip().split("\t")
            qseqid = line[0]
            rseqid = line[1]

            if qseqid == rseqid:
                # for self-blast searches to find paralogs, skip self-hits
                continue
            try:
                bitscore = float(line[-1])
            except:
                raise RuntimeError(f"{line} \nfrom {blast_infile_path}\ncould not be parsed\n")
            if qseqid not in best_hits_dict:
                best_hits_dict[qseqid] = BestHit(qseqid=qseqid, rseqid=rseqid, bitscore=bitscore)
            else:
                best_hits_dict[qseqid].update_besthit(new_qseqid = qseqid, new_rseqid = rseqid, new_bitscore = bitscore)
    return best_hits_dict  


def get_BRHs(besthits_infile1, besthits_infile2, annotation1, annotation2, x_list1, x_list2, y_list1 = [], y_list2 = [], species1 = "", species2 = "", outfile_path:str = "", contig_names=False):
    """
    get a dictionary with {species1_ID : species2_ID} of all best reciprocal hits
    """
    out_dict = {}
    if type(annotation1) == str and type(annotation2)== str:
        ## use annotation filenames as species headers
        species1 = annotation1.split("/")[-1].split(".")[0]
        species2 = annotation2.split("/")[-1].split(".")[0]
        ## read annotations from filepaths
        print(f"reading annotations of species 1 and 2...")
        try:
            annotation1 = gff.parse_gff3_general(annotation1, verbose=False)
        except:
            print(f"-- > parse {annotation1} as gtf")
            annotation1 = gff.parse_gff3_general(annotation1, verbose=True, gtf=True)
        try:
            annotation2 = gff.parse_gff3_general(annotation2, verbose=False)
        except:
            print(f"-- > parse {annotation2} as gtf")
            annotation2 = gff.parse_gff3_general(annotation2, verbose=True, gtf=True)

    elif type(annotation1) == dict and type(annotation2) == dict:
        ## use species names from parameters and assume that the annotations are already read in
        if species1 == "" or species2 == "":
            raise RuntimeError(f"if you include parsed annotations you have to give species names in the function parameters. You have given:\n species1 = '{species1}'\n species2 = '{species2}'")
        pass

    if contig_names:
        header = f"{species1}\tchromosome\tcontig_ID\t{species2}\tchromosome\tcontig_ID\n"
    else:
        header = f"{species1}\tchromosome\t{species2}\tchromosome\n"

    for species1_id, species1_besthit in besthits_infile1.items():
        species2_id = species1_besthit.rseqid
        try:
            species2_besthit = besthits_infile2[species2_id]
        except:
            continue
        if species1_id == species2_besthit.rseqid:
            out_dict[species1_id] = species2_id
    if outfile_path != "":
        if len(out_dict)>0:
            with open(outfile_path, "w") as outfile:
                outfile.write(header)
                for species1_id, species2_id in out_dict.items():
                    ## remove the "_1" suffix so that they can be found in the annotation
                    if species1_id[-2:] == "_1":
                        species1_id = species1_id[:-2]
                    if species2_id[-2:] == "_1":
                        species2_id = species2_id[:-2]
                    contig1 = "A"
                    contig2 = "A"
                    try:
                        ID_contig1 = annotation1[species1_id].contig
                    except:
                        raise RuntimeError(f"{species1_id} in {species1} not found in the annotation")
                    try:
                        ID_contig2 = annotation2[species2_id].contig
                    except:
                        raise RuntimeError(f"{species2_id} in {species2} not found in the annotation")
                    if ID_contig1 in x_list1:
                        contig1 = "X"
                    elif ID_contig1 in y_list1:
                        contig1 = "Y"
                    if ID_contig2 in x_list2:
                        contig2 = "X"
                    elif ID_contig2 in y_list2:
                        contig2 = "Y"
                    if contig_names:
                        outfile.write(f"{species1_id}\t{contig1}\t{ID_contig1}\t{species2_id}\t{contig2}\t{ID_contig2}\n")
                    else:
                        outfile.write(f"{species1_id}\t{contig1}\t{species2_id}\t{contig2}\n")
        print(f"outfile written to: {outfile_path}")
    return out_dict


if __name__ == "__main__":
    

    args=parse_args()
    blast_infile_path1 = args.blast1
    blast_infile_path2 = args.blast2
    outfile = args.outfile
    annotation_path1 = args.annotation1
    annotation_path2 = args.annotation2
    verbose = args.verbose
    ###
    # the C. magnifica contigs contain a comma because God is trying to test me. don't split in this case
    ###
    if not args.X_contigs1:
        sex_chr_contigs = sex_chromosome_names()
        species1_listname = blast_infile_path1.split("/")[-1]
        species1 = gff.split_at_second_occurrence(species1_listname)
        species2_listname = blast_infile_path2.split("/")[-1]
        species2 = gff.split_at_second_occurrence(species2_listname)
        X_list1 = sex_chr_contigs[species1]["X"]
        X_list2 = sex_chr_contigs[species2]["X"]
        Y_list1 = sex_chr_contigs[species1]["Y"]
        Y_list2 = sex_chr_contigs[species2]["Y"]
        # autosomes = autosomes_lists()
        # A_list1 = autosomes[species1]
        # A_list2 = autosomes[species2]
    else:
        X_list1 = args.X_contigs1.strip().split(",")
        X_list2 = args.X_contigs2.strip().split(",")
        if args.Y_contigs1 and args.Y_contigs2:
            Y_list1 = args.Y_contigs1.strip().split(",")
            Y_list2 = args.Y_contigs2.strip().split(",")
        else:
            Y_list1 = []
            Y_list2 = []

    besthits_infile1 = read_best_hits(blast_infile_path1)
    besthits_infile2 = read_best_hits(blast_infile_path2)
    # print(besthits_infile1["rna-AOBTE_LOCUS3-2_1"])

    BRH_dict = get_BRHs(besthits_infile1, besthits_infile2, annotation1=annotation_path1, annotation2=annotation_path2, 
        x_list1=X_list1, x_list2=X_list2, 
        y_list1=Y_list1, y_list2=Y_list2, 
        contig_names=True,
        outfile_path=outfile)

