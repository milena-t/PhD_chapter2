"""
filter all the cross intersections to only keep the ones that are the best in each pair
"""

import sys,random
import miniprot_stats_comparison as minialn

if __name__=="__main__":
    overlap_list = sys.argv[1]
    miniprot_annot_path = sys.argv[2]
    outlist = sys.argv[3]

    miniprot_dict = minialn.miniprot_parse_alignment(miniprot_annot_path, queryIDs=False)
    # print(f"------")
    # print(miniprot_dict["MP027752"])
    # print(miniprot_dict["MP028525"])
    # print(f"------")

    overlap_other=0
    overlap_list_worse = []
    overlap_matches = {}
    with open(overlap_list, "r") as overlap_lines, open(outlist, "w") as outlist_file:
        for line in overlap_lines.readlines():
            id1, id2 = line.strip().split()
            if id1==id2:
                continue
            else:
                overlap_other+=1
                overlap_matches.setdefault(id1, set()).add(id2)
                overlap_matches.setdefault(id2, set()).add(id1)

        # sort ids by their "scores" (aligment seq ident, alignment length)
        scores = {} 
        for id_, aln in miniprot_dict.items():
            scores[id_] = (aln.identity, aln.length)
        ids = list(scores)
        random.shuffle(ids) # randomize in case of ties
        order = sorted(ids, key=lambda x: scores[x], reverse=True)

        # loop through all IDs starting from the best, if an ID is hit for the first time, 
        # it is clearly the top of its node according to the score sorting, so it is kept, 
        # and all other matches (dict value list members) are removed.
        # if an id is hit that was removed (is an inferior member of an overlap cluster)
        # it is ignored since it was already handled
        removed, kept = set(), []
        for id_ in order:
            if id_ in removed:
                continue
            kept.append(id_)
            removed.update(overlap_matches.get(id_, ()))

        outlist_file.write("\n".join(kept))

    print(f"{len(ids)} overlap with non-self mapping queries (reciprocal)")
    print(f"{len(kept)} clusters, {len(removed)} IDs removed, ({len(kept)} + {len(removed)} = {len(kept)+len(removed)})")
    


