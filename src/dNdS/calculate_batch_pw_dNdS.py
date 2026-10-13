"""
Make batches of dNdS to calculate at once. splits it into several jobs of orthologs
"""

import os
import sys


def make_nested_lists(dir_path, include_list = []):
    """
    make dictionary of by-species nested lists for the pairwise fasta files. if include_list is specified then only use 
    pairs that have members in that list
    """
    subdirs = [f"{dir_path}{f}/" for f in os.listdir(dir_path) if os.path.isdir(os.path.join(dir_path, f))]
    if include_list != []:
        print(f"include only species: {include_list}")

    fastas = []
    for pairdir in subdirs:
        print(f"making nested fasta lists in: {pairdir}")
        species_present = True
        if include_list != []:
            species_present = False
            for species_incl in include_list:
                print(f"\t - testing {species_incl} in {pairdir}")
                if species_incl in pairdir:
                    species_present = True
                    break
        if species_present:
            fastas.extend([f"{pairdir}{f}" for f in os.listdir(pairdir) if os.path.isfile(os.path.join(pairdir, f))])
    
    nested_dict = {}
    print(f"{len(fastas)} orthologs: {fastas[:3]} ... ")
        
    for fasta_name in fastas:
        fasta_components = fasta_name.split("/")[-1].split("_")
        species1 = f"{fasta_components[0]}_{fasta_components[1]}"
        species2 = f"{fasta_components[2]}_{fasta_components[3]}"
        if species1 in nested_dict:
            if species2 in nested_dict[species1]:
                nested_dict[species1][species2].append(fasta_name)
            else:
                nested_dict[species1][species2] = [fasta_name]
        else:
            nested_dict[species1] = {}
            nested_dict[species1][species2] = [fasta_name]

    return nested_dict


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    out_list = []
    for i in range(0, len(lst), n):
        out_list.append(lst[i:i + n])
    return out_list



if __name__ == "__main__":

    # list files in dir
    datadir = f"/Users/miltr339/work/pairwise_blast_chapter_2_3/brh_tables/"
    bash_dir = f"/Users/miltr339/work/PhD_code/PhD_chapter3/bash"
    dNdS_exec = f"bash"
    pelle = False

    if not os.path.isdir(datadir):
        
        # revisions
        datadir = f"/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/dNdS_diorhabda_X/"

        bash_dir = f"/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/PhD_chapter2/bash"
        dNdS_exec = f"sbatch"
        pelle = True
    
    #############################
    ###############
    #######
    
    if len(sys.argv)==2:
        chr_type = sys.argv[1]
    else:
        chr_type = "A"
    
    #######
    
    # analysis = "LRT" ## codeml site model M1a and M2a with likelihood ratio test
    # analysis = "dNdS" ## codeml branch model
    
    #### revision
    # analysis = "M7-M8_LRT"
    # analysis = "LRT"
    analysis = "dNdS" ## codeml branch model
    
    #######
    ###############
    #############################
    ## revision: 

    script_path = f"{bash_dir}/run_batch_{analysis}_chrysomelidae.sh"
    if analysis == "M7-M8_LRT":
        chr_type = "A"
        outdir = f"{datadir}site_model_beta_res_{chr_type}/"
    elif analysis == "LRT":
        chr_type = "A"
        outdir = f"{datadir}site_model_res_{chr_type}/"
    
    #### chrysomelidae
    elif analysis == "dNdS":
        chr_type = "A"
        outdir = f"{datadir}branch_model_results_{chr_type}/"

    
    if True:
        ########
        ## for testing purposes only do a few files of each pair
        ## if 0 then it takes all files, for a normal run
        num_files = 0
        ########
        os.chdir(outdir)
        print(f"ANALYZING : {analysis}")
                
        ### split into several batches because it doesn't run with too many fasta files
        print(datadir)
        fastadir=f"{datadir}multifasta_diorhabdaX_bruchini{chr_type}/"
        print(fastadir)
        fasta_overall_list = [f"{fastadir}{f}" for f in os.listdir(fastadir)  if ".fasta" in f]
        
        fasta_unique_list = sorted(list(set(fasta_overall_list)))
        print(f"\n{len(fasta_unique_list)} unique fasta files: {fasta_unique_list[:2]}...\n")
        
        batch_size = 100
        i = 0 # counter

        sep_lists = chunks(fasta_unique_list, batch_size)
        print(f"run {len(sep_lists)} jobs ... \n")
        for fasta_list_sublist in sep_lists:

            # re-do specific batches, 4,5 or 20,25 where there is something weird
            if i not in [2]:
                print("------------------------------")
                print(f"skip : {i}")
                i +=1
                continue

            # fasta_list_sublist.reverse() 
            # some batches contain orthologs where the loop does not continue, but the job also does not terminate
            # So I re-run those from the other side (reversed list) to get all orthologs in that batch. 
            # remember to rename the first output files to keep the results from the first half of the orthologs in this batch!
            
            print("------------------------------")
            print(f"\t * {len(fasta_list_sublist)} fasta list: {fasta_list_sublist[0][51:]}...{fasta_list_sublist[-1][51:]}")
            fasta_string = " ".join(fasta_list_sublist)
            if pelle:
                jobname = f"{i}_{analysis}_{chr_type}-linked"
                dNdS_command = f"{dNdS_exec} -J {jobname} -o {jobname}.out {script_path} {fasta_string}"
            else:
                dNdS_command = f"{dNdS_exec} {script_path} {fasta_string}"

            if True:
                os.system(dNdS_command)
                os.system("sleep .5") # wait a little bit after each command so that the job manager can keep up
            comm_spl=dNdS_command.split(" ")
            print(len(comm_spl))
            dNdS_command_print="\n".join(comm_spl)
            print(f"\t{dNdS_command_print[:500]}...{dNdS_command_print[-50:]}\n")
            i +=1
        print(f"FILES IN : {outdir}")

