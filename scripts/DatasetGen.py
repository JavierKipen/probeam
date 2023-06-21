# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 14:35:57 2023

@author: JK-WORK
"""

from cleave_proteins import cleave_proteins
from dye_seqs_from_peptides import dye_seqs_from_peptides
import os, shutil
import subprocess


def erase_contents(folder): #Erase contents of folder
    for filename in os.listdir(folder): ## https://stackoverflow.com/questions/185936/how-to-delete-the-contents-of-a-folder
        file_path = os.path.join(folder, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print('Failed to delete %s. Reason: %s' % (file_path, e))

n_proteins=[10,20,50,100,200,500,1000,2000,5000,10000,20000];
#n_proteins=[2000,5000,20000];

path_datasets_figure="../Own/HMM_modif/Datasets/ForPaper/"
path_datasets_figure_common_files= path_datasets_figure +"CommonFiles/"
prot_fasta=path_datasets_figure_common_files+"UP000005640_9606.fasta"
seq_params_path= path_datasets_figure_common_files + "seq-params.json"

for n in n_proteins:
    protein_folder=path_datasets_figure+str(n)+"Prot"
    if not os.path.isdir(protein_folder): ##If folder already Exists
        os.mkdir(protein_folder)

        #erase_contents(protein_folder)
    #else:
    
    protein_folder += "/"
    
    cleave_proteins(prot_fasta,
                protein_folder+"peptides.csv",
                "trypsin",
                n=n)
    dye_seqs_path=protein_folder+"dye-seqs.tsv";
    dye_seqs_from_peptides(protein_folder+"peptides.csv",
                       ['DE','C','Y'],
                       dye_seqs_path)
    dye_tracks_path=protein_folder+"dye-tracks.tsv";
    radiometries_path=protein_folder+"radiometries.tsv";
    true_ids_path=protein_folder+"true-ids.tsv";
    cmd_gen_dye_tracks = "./../cc_code/bin/release/whatprot simulate dt -t 10 -g 1000 -P " + seq_params_path + " -S " + dye_seqs_path + " -T " +dye_tracks_path;
    
    cmd_sim_rad = "./../cc_code/bin/release/whatprot simulate rad -t 10 -g 100000 -P " + seq_params_path + " -S " + dye_seqs_path + " -R " + radiometries_path + " -Y "+ true_ids_path 
    
    subprocess.run(cmd_sim_rad, shell=True, check=True)
    subprocess.run(cmd_gen_dye_tracks, shell=True, check=True)
    #subprocess.run(cmd_sim_rad)
