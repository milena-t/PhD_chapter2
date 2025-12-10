"""
Get the transcript IDs of best reciprocal hits (BRH) of two species. the input is the outfmt 6 blast results in either direction
useage: python3 get_blast_BRH.py species1_species2.blast species2_species1.blast
and the output is a tsv file with two columns, one is species1 and one is species2, and each row contains a set of gene IDs that are each others BRH
I identify the best hit by the highest bit score
"""
import sys

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
            bitscore = float(line[-1])
            if qseqid not in best_hits_dict:
                best_hits_dict[qseqid] = BestHit(qseqid=qseqid, rseqid=rseqid, bitscore=bitscore)
            else:
                best_hits_dict[qseqid].update_besthit(new_qseqid = qseqid, new_rseqid = rseqid, new_bitscore = bitscore)
    return best_hits_dict  


def get_BRHs(besthits_infile1, besthits_infile2, outfile_path:str = ""):
    """
    get a dictionary with {species1_ID : species2_ID} of all best reciprocal hits
    """
    out_dict = {}
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
                for species1_id, species2_id in out_dict.items():
                    outfile.write(f"{species1_id}\t{species2_id}\n")
        print(f"outfile written to: {outfile_path}")
    return out_dict


if __name__ == "__main__":
    
    blast_infile_path1 = sys.argv[1]
    blast_infile_path2 = sys.argv[2]
    if len(sys.argv)==3:
        outfile = f"{blast_infile_path1}_BRH.tsv"
    elif len(sys.argv) == 4:
        outfile_name = sys.argv[3]

    besthits_infile1 = read_best_hits(blast_infile_path1)
    besthits_infile2 = read_best_hits(blast_infile_path2)
    print(besthits_infile1["rna-AOBTE_LOCUS3-2_1"])

    BRH_dict = get_BRHs(besthits_infile1, besthits_infile2, outfile_path=outfile)