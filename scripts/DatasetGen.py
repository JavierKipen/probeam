# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 14:35:57 2023

@author: JK-WORK
"""

import os, shutil,sys
import subprocess

import sys, os
path_to_whatprot_python=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+ '/ext/whatprot/python';
sys.path.append(path_to_whatprot_python)
from cleave_proteins import cleave_proteins
from dye_seqs_from_peptides import dye_seqs_from_peptides


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

n_proteins=[20,50,100,200,500,1000,2000,5000,10000,20000];


path_datasets=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+"/data/"
path_datasets_norm=path_datasets+"NormDatasets/"
path_datasets_long=path_datasets+"LongDatasets/"
path_datasets_common_files= path_datasets +"common/"
prot_fasta=path_datasets_common_files+"UP000005640_9606.fasta"
seq_params_path= path_datasets_common_files + "seq-params.json"

n_reads_norm=100000;
n_reads_long=1000000;


if not os.path.isdir(path_datasets_norm):
    os.mkdir(path_datasets_norm)
if not os.path.isdir(path_datasets_long):
    os.mkdir(path_datasets_long)

##Generate norm datasets
for n in n_proteins:
    protein_folder=path_datasets_norm+str(n)+"Prot/"
    if not os.path.isdir(protein_folder): ##If folder doesnt exist
        os.mkdir(protein_folder)
    
    cleave_proteins(prot_fasta,protein_folder+"peptides.csv","trypsin",n=n)
    dye_seqs_path=protein_folder+"dye-seqs.tsv";
    dye_seqs_from_peptides(protein_folder+"peptides.csv",
                       ['DE','C','Y'],
                       dye_seqs_path)
    dye_tracks_path=protein_folder+"dye-tracks.tsv";
    radiometries_path=protein_folder+"radiometries.tsv";
    true_ids_path=protein_folder+"true-ids.tsv";
    cmd_gen_dye_tracks = "./bin/whatprot simulate dt -t 10 -g 1000 -P " + seq_params_path + " -S " + dye_seqs_path + " -T " +dye_tracks_path;
    
    cmd_sim_rad = "./bin/whatprot simulate rad -t 10 -g "+str(n_reads_norm)+" -P " + seq_params_path + " -S " + dye_seqs_path + " -R " + radiometries_path + " -Y "+ true_ids_path 
    
    subprocess.run(cmd_sim_rad, shell=True, check=True)
    subprocess.run(cmd_gen_dye_tracks, shell=True, check=True)
    #subprocess.run(cmd_sim_rad)

#Generate long datasets:
n_proteins_long=1000
protein_folder=path_datasets_long+str(n_proteins_long)+"Prot/"
if not os.path.isdir(protein_folder): ##If folder doesnt exist
    os.mkdir(protein_folder)

cleave_proteins(prot_fasta,protein_folder+"peptides.csv","trypsin",n=n_proteins_long)
dye_seqs_path=protein_folder+"dye-seqs.tsv";
dye_seqs_from_peptides(protein_folder+"peptides.csv",
                       ['DE','C','Y'],
                       dye_seqs_path)
dye_tracks_path=protein_folder+"dye-tracks.tsv";
radiometries_path=protein_folder+"radiometries.tsv";
true_ids_path=protein_folder+"true-ids.tsv";
cmd_gen_dye_tracks = "./bin/whatprot simulate dt -t 10 -g 1000 -P " + seq_params_path + " -S " + dye_seqs_path + " -T " +dye_tracks_path;

cmd_sim_rad = "./bin/whatprot simulate rad -t 10 -g "+str(n_reads_long)+" -P " + seq_params_path + " -S " + dye_seqs_path + " -R " + radiometries_path + " -Y "+ true_ids_path 

subprocess.run(cmd_sim_rad, shell=True, check=True)
subprocess.run(cmd_gen_dye_tracks, shell=True, check=True)
