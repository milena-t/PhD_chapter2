"""
Loop circos through all the karyotypes and species automatically
Needs to modify circos.conf for input files and image.conf for output file
"""

import os 
import subprocess

def get_karyotype_files_cmac_populations(username="miltr339"):
    filesdir =f"/Users/{username}/work/PhD_code/PhD_chapter2/data/circos/Cmac_populations"
    outdir = {
        "China" : {
            "Lome" : f"{filesdir}/China_Lome_circos_karyotype.txt",
            "China" : f"{filesdir}/China_China_circos_karyotype.txt",
        },
        "Lome" : {
            "China" : f"{filesdir}/Lome_China_circos_karyotype.txt",
        }
    }
    return outdir


def get_karyotype_files(username="miltr339"):
    filesdir =f"/Users/{username}/work/PhD_code/PhD_chapter2/data/circos/"
    outdir = {
        "A_obtectus" : {
            "A_obtectus" : f"{filesdir}A_obtectus_A_obtectus_circos_karyotype.txt",
            "B_siliquastri" : f"{filesdir}A_obtectus_B_siliquastri_circos_karyotype.txt",
            "B_varius" : f"{filesdir}A_obtectus_B_varius_circos_karyotype.txt",
            "C_chinensis" : f"{filesdir}A_obtectus_C_chinensis_circos_karyotype.txt",
            "C_maculatus" : f"{filesdir}A_obtectus_C_maculatus_circos_karyotype.txt",
            "D_carinulata" : f"{filesdir}A_obtectus_D_carinulata_circos_karyotype.txt",
            "D_sublineata" : f"{filesdir}A_obtectus_D_sublineata_circos_karyotype.txt",
        },
        "B_siliquastri" : {
            "A_obtectus" : f"{filesdir}B_siliquastri_A_obtectus_circos_karyotype.txt",
            "B_siliquastri" : f"{filesdir}B_siliquastri_B_siliquastri_circos_karyotype.txt",
            "B_varius" : f"{filesdir}B_siliquastri_B_varius_circos_karyotype.txt",
            "C_chinensis" : f"{filesdir}B_siliquastri_C_chinensis_circos_karyotype.txt",
            "C_maculatus" : f"{filesdir}B_siliquastri_C_maculatus_circos_karyotype.txt",
            "D_carinulata" : f"{filesdir}B_siliquastri_D_carinulata_circos_karyotype.txt",
            "D_sublineata" : f"{filesdir}B_siliquastri_D_sublineata_circos_karyotype.txt",
        },
        "B_varius" : {
            "A_obtectus" : f"{filesdir}B_varius_A_obtectus_circos_karyotype.txt",
            "B_siliquastri" : f"{filesdir}B_varius_B_siliquastri_circos_karyotype.txt",
            "B_varius" : f"{filesdir}B_varius_B_varius_circos_karyotype.txt",
            "C_chinensis" : f"{filesdir}B_varius_C_chinensis_circos_karyotype.txt",
            "C_maculatus" : f"{filesdir}B_varius_C_maculatus_circos_karyotype.txt",
            "D_carinulata" : f"{filesdir}B_varius_D_carinulata_circos_karyotype.txt",
            "D_sublineata" : f"{filesdir}B_varius_D_sublineata_circos_karyotype.txt",
        },
        "C_chinensis" : {
            "A_obtectus" : f"{filesdir}C_chinensis_A_obtectus_circos_karyotype.txt",
            "B_siliquastri" : f"{filesdir}C_chinensis_B_siliquastri_circos_karyotype.txt",
            "B_varius" : f"{filesdir}C_chinensis_B_varius_circos_karyotype.txt",
            "C_chinensis" : f"{filesdir}C_chinensis_C_chinensis_circos_karyotype.txt",
            "C_maculatus" : f"{filesdir}C_chinensis_C_maculatus_circos_karyotype.txt",
            "D_carinulata" : f"{filesdir}C_chinensis_D_carinulata_circos_karyotype.txt",
            "D_sublineata" : f"{filesdir}C_chinensis_D_sublineata_circos_karyotype.txt",
        },
        "C_maculatus" : {
            "A_obtectus" : f"{filesdir}C_maculatus_A_obtectus_circos_karyotype.txt",
            "B_siliquastri" : f"{filesdir}C_maculatus_B_siliquastri_circos_karyotype.txt",
            "B_varius" : f"{filesdir}C_maculatus_B_varius_circos_karyotype.txt",
            "C_chinensis" : f"{filesdir}C_maculatus_C_chinensis_circos_karyotype.txt",
            "C_maculatus" : f"{filesdir}C_maculatus_C_maculatus_circos_karyotype.txt",
            "D_carinulata" : f"{filesdir}C_maculatus_D_carinulata_circos_karyotype.txt",
            "D_sublineata" : f"{filesdir}C_maculatus_D_sublineata_circos_karyotype.txt",
        },
        "D_carinulata" : {
            "A_obtectus" : f"{filesdir}D_carinulata_A_obtectus_circos_karyotype.txt",
            "B_siliquastri" : f"{filesdir}D_carinulata_B_siliquastri_circos_karyotype.txt",
            "B_varius" : f"{filesdir}D_carinulata_B_varius_circos_karyotype.txt",
            "C_chinensis" : f"{filesdir}D_carinulata_C_chinensis_circos_karyotype.txt",
            "C_maculatus" : f"{filesdir}D_carinulata_C_maculatus_circos_karyotype.txt",
            "D_carinulata" : f"{filesdir}D_carinulata_D_carinulata_circos_karyotype.txt",
            "D_sublineata" : f"{filesdir}D_carinulata_D_sublineata_circos_karyotype.txt",
        },
        "D_sublineata" : {
            "A_obtectus" : f"{filesdir}D_sublineata_A_obtectus_circos_karyotype.txt",
            "B_siliquastri" : f"{filesdir}D_sublineata_B_siliquastri_circos_karyotype.txt",
            "B_varius" : f"{filesdir}D_sublineata_B_varius_circos_karyotype.txt",
            "C_chinensis" : f"{filesdir}D_sublineata_C_chinensis_circos_karyotype.txt",
            "C_maculatus" : f"{filesdir}D_sublineata_C_maculatus_circos_karyotype.txt",
            "D_carinulata" : f"{filesdir}D_sublineata_D_carinulata_circos_karyotype.txt",
            "D_sublineata" : f"{filesdir}D_sublineata_D_sublineata_circos_karyotype.txt",
        },
    }
    return outdir

def modify_circos_conf(circos_conf, config_files_dict, verbose=True):
    """
    modify the config file
    """
    modified_lines = []        
    with open(circos_conf, "r") as codeml:
        lines = codeml.readlines()
        for line in lines: # go through all lines
            for key, value in config_files_dict.items():
                if key in line: # check if to-modify variable is in line
                    line = f"{key}{value}\n" # make new line and overwrite the old one
                    if verbose:
                        print("\t"+line.split("\n")[0]) # remove the newline character for printing so it looks nicer
            modified_lines.append(line)
    with open(circos_conf, "w") as config:
        config.writelines(modified_lines)

    if verbose:
        print(f"\t\t--> done modifying {circos_conf}")
        print()


if __name__ == "__main__":
    username = "milena"
    # karyotype_dict = get_karyotype_files(username=username)
    karyotype_dict = get_karyotype_files_cmac_populations(username=username)
    species_list = list(karyotype_dict.keys())

    circos_infiles_dir = f"/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/PhD_chapter2/data/circos"
    # circos_infiles_dir = f"/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/PhD_chapter2/data/circos/Cmac_populations"
    outfiles_dir = f"/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/circos/plots"


    os.chdir(outfiles_dir)
    for species1 in species_list:
        for species2 in species_list:
            print(f" ===================== {species1} vs. {species2} =====================")
            infiles_dict = {
                "karyotype = " : karyotype_dict[species1][species2],
                # "file          = " : f"{circos_infiles_dir}/circos_links_{species1}_vs_{species2}.txt",
                "file          = " : f"{circos_infiles_dir}/Cmac_populations/circos_links_{species1}_vs_{species2}.txt",
            }
            modify_circos_conf(circos_conf=f"{circos_infiles_dir}/conf_files/circos.conf", config_files_dict=infiles_dict)
            outfile = {
                "dir   = " : f".", # f"{outfiles_dir}",
                "file  = " : f"{species1}_{species2}_circos.png",
            }
            modify_circos_conf(circos_conf=f"{circos_infiles_dir}/conf_files/image.conf", config_files_dict=outfile)

            ## TODO
            conf_path = f"{circos_infiles_dir}/conf_files/circos.conf"
            inlist = ["circos", "-conf", conf_path]
            subprocess.run(inlist, check=True)#,  stdout=subprocess.DEVNULL)
            print(f" ".join(inlist))

            for chr in ["X", "Y"]:
                infiles_dict = {
                    "karyotype = " : karyotype_dict[species1][species2],
                    # "file          = " : f"{circos_infiles_dir}/circos_links_{species1}_vs_{species2}_{chr}.txt",
                    "file          = " : f"{circos_infiles_dir}/Cmac_populations/circos_links_{species1}_vs_{species2}_{chr}.txt",
                }
                modify_circos_conf(circos_conf=f"{circos_infiles_dir}/conf_files/circos.conf", config_files_dict=infiles_dict)
                outfile = {
                    "dir   = " : f".", # f"{outfiles_dir}",
                    "file  = " : f"{species1}_{species2}_circos_{chr}.png",
                }
                modify_circos_conf(circos_conf=f"{circos_infiles_dir}/conf_files/image.conf", config_files_dict=outfile)

                ## TODO
                conf_path =f"{circos_infiles_dir}/conf_files/circos.conf"
                inlist = ["circos", "-conf", conf_path]
                subprocess.run(inlist, check=True)#,  stdout=subprocess.DEVNULL)
                print(f" ".join(inlist))

            # break
        # break

