"""
Loop circos through all the karyotypes and species automatically
Needs to modify circos.conf for input files and image.conf for output file
"""

import os 
import subprocess

def get_karyotype_files_cmac_populations(username="miltr339"):
    filesdir =f"/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/PhD_chapter2/data/circos/Cmac_populations"
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
    filesdir =f"/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/PhD_chapter2/data/circos/"
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


    karyotype_dict = get_karyotype_files(username=username)
    circos_infiles_dir = f"/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/PhD_chapter2/data/circos"

    # karyotype_dict = get_karyotype_files_cmac_populations(username=username)
    # circos_infiles_dir = f"/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/PhD_chapter2/data/circos/Cmac_populations"

    species_list = list(karyotype_dict.keys())
    outfiles_dir = f"/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/circos/plots"
    # rsync -azP "milenatr@pelle.uppmax.uu.se:/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/circos/plots/*png" /Users/miltr339/work/PhD_code/PhD_chapter2/data/circos/plots

    nucleotide_blast = False
    
    
    os.chdir(outfiles_dir)
    for species1 in species_list:
        for species2 in species_list:

            if species1=="Lome" and species2=="Lome" :
                continue

            if species1 != species2:
                continue

            print(f" ===================== {species1} vs. {species2} =====================")
            if "Lome" in species_list:
                infiles_dict = {
                    "karyotype = " : karyotype_dict[species1][species2],
                    "file          = " : f"{circos_infiles_dir}/Cmac_populations/circos_links_{species1}_vs_{species2}.txt",
                }
            else:
                infiles_dict = {
                    "karyotype = " : karyotype_dict[species1][species2],
                    "file          = " : f"{circos_infiles_dir}/circos_links_{species1}_vs_{species2}.txt",
                }
            outfile = {
                "dir   = " : f".", # f"{outfiles_dir}",
                "file  = " : f"{species1}_{species2}_circos.png",
            }
                
            if nucleotide_blast and species1==species2:
                # this cannot be plotted with nucleotide blast since it hits the circos max links limit
                pass
            else:

                print(f" -------------------- all chromosomes --------------------")
                modify_circos_conf(circos_conf=f"{circos_infiles_dir}/conf_files/circos.conf", config_files_dict=infiles_dict)
                modify_circos_conf(circos_conf=f"{circos_infiles_dir}/conf_files/image.conf", config_files_dict=outfile)

                conf_path = f"{circos_infiles_dir}/conf_files/circos.conf"
                inlist = ["circos", "-conf", conf_path]
                subprocess.run(inlist, check=True,  stdout=subprocess.DEVNULL)
                print(f" ".join(inlist))

            for chr in ["X", "Y"]:
                print(f" -------------------- {chr}-chromosome --------------------")
                if "Lome" in species_list:
                    infiles_dict = {
                        "karyotype = " : karyotype_dict[species1][species2],
                        "file          = " : f"{circos_infiles_dir}/Cmac_populations/circos_links_{species1}_vs_{species2}_{chr}.txt",
                    }
                else:
                    infiles_dict = {
                        "karyotype = " : karyotype_dict[species1][species2],
                        "file          = " : f"{circos_infiles_dir}/circos_links_{species1}_vs_{species2}_{chr}.txt",
                    }
                outfile = {
                    "dir   = " : f".", # f"{outfiles_dir}",
                    "file  = " : f"{species1}_{species2}_circos_{chr}.png",
                }

                if nucleotide_blast and species1==species2:
                    print(f" ===================== {species1} nucleotide blast =====================")
                    infiles_dict = {
                        "karyotype = " : karyotype_dict[species1][species2],
                        "file          = " : f"{circos_infiles_dir}/circos_links_{species1}_nucleotide_blast_{chr}.txt",
                    }
                    outfile = {
                        "dir   = " : f".", # f"{outfiles_dir}",
                        "file  = " : f"{species1}_nucleotide_blast_circos_{chr}.png",
                    }
                else:

                    # add link density histogram
                    try:
                        conf_path =f"{circos_infiles_dir}/conf_files/circos_hist.conf"
                        link_density_hist_exe=f"/sw/apps/circos/0.69-9/rackham/circos-tools-0.23/tools/binlinks/bin/binlinks"
                
                        linksfile = infiles_dict["file          = "]
                        # hist_outfile_ = linksfile.replace(".txt", "_hist_vals.txt").split("/")[-1]
                        hist_outfile_ = f"circos_links_{species1}_vs_{species2}_{chr}_hist_vals.txt" 
                        hist_outfile = f"/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/circos/histograms/{hist_outfile_}"
                        command = [link_density_hist_exe, "-links", linksfile , ">", hist_outfile]
                        print(" ".join(command))
                        print(f"--> make histogram file:")
                        # with open(hist_outfile, "w") as f:
                        #     subprocess.run(command, check=True, stdout=f)
                        subprocess.run(command, check=True, shell=True)
                        if os.path.isfile(hist_outfile):
                            print(f"successfully created histogram file: {hist_outfile}")
                        else:
                            raise RuntimeError(f"**!!  histogram file was not successfully created with command:\n{command}")
                        infiles_dict["file      = "] = hist_outfile
                    except:
                        conf_path =f"{circos_infiles_dir}/conf_files/circos.conf"
                        print(f"no link density histogram could be generated for {linksfile}")
                        pass

                    modify_circos_conf(circos_conf=conf_path, config_files_dict=infiles_dict)
                    modify_circos_conf(circos_conf=f"{circos_infiles_dir}/conf_files/image.conf", config_files_dict=outfile)

                    inlist = ["circos", "-conf", conf_path]
                    try:
                        subprocess.run(inlist, check=True,  stdout=subprocess.DEVNULL)
                        print(f" ".join(inlist))
                    except:
                        continue

            # break
        # break

