# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 17:10:28 2023

@author: JK-WORK
"""

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

def get_proc_time(outputHybrid):
    output=str(outputHybrid);
    idx_t_start=output.find("Finished classification")+25;
    substr=output[idx_t_start:];
    computing_time=float(substr[:substr.find(" ")])
    return computing_time;

def get_proc_time_beam(outputBeam):
    output=str(outputBeam);
    idx_t_start=output.find(": ")+2;
    substr=output[idx_t_start:];
    computing_time=float(substr[:substr.find(" ")])
    return computing_time;

n_proteins=[10,20,50,100,200,500,1000,2000,5000,10000,20000];
#n_proteins=[20000];

path_datasets_figure="../Own/HMM_modif/Datasets/ForPaper/"
path_datasets_figure_common_files= path_datasets_figure +"CommonFiles/"
prot_fasta=path_datasets_figure_common_files+"UP000005640_9606.fasta"
seq_params_path= path_datasets_figure_common_files + "seq-params.json"

full_path_to_ds="~/Fluoroseq/Own/HMM_modif/Datasets/ForPaper/"

cmd_set_threads="export OMP_NUM_THREADS=1";
subprocess.run(cmd_set_threads, shell=True, check=True) ##Sets to only one thread running

for n in n_proteins:
    print("Number of proteins: "+ str(n))
    protein_folder=path_datasets_figure+str(n)+"Prot/"
    dye_seqs_path=protein_folder+"dye-seqs.tsv";
    dye_tracks_path=protein_folder+"dye-tracks.tsv";
    radiometries_path=protein_folder+"radiometries.tsv";
    true_ids_path=protein_folder+"true-ids.tsv";
    predictions_hybrid_path = protein_folder+"predictionsHybrid.csv"
    
    cmd_classify_hybrid = "./../cc_code/bin/release/whatprot classify hybrid -k 10000 -s 0.5 -H 1000 -p 5 -P " + seq_params_path + " -S " + dye_seqs_path + " -T " +dye_tracks_path + " -R " + radiometries_path + " -Y " + predictions_hybrid_path;
    
    cmd_classify_beam = "./../Own/HMM_modif/LowLevel/Rewriting/BeamDec " + full_path_to_ds + str(n)+"Prot/";
    
    #output=subprocess.run(cmd_classify_hybrid, shell=True, check=True)
    outputHybrid=subprocess.check_output("export OMP_NUM_THREADS=1\n"+cmd_classify_hybrid, shell=True)
    computing_time_hybrid=get_proc_time(outputHybrid)
    hybrid_msg="Hybrid computing total time (seconds): "+ str(computing_time_hybrid);
    print(str(outputHybrid))
    outputBeam=subprocess.check_output(cmd_classify_beam, shell=True)
    computing_time_beam=get_proc_time_beam(outputBeam)
    beam_msg="Beam computing time per read (miroseconds): "+ str(computing_time_beam)
    print(str(outputBeam))
    with open(protein_folder+"timing-results.txt", 'w') as f:
        f.write(hybrid_msg)
        f.write('\n')
        f.write(beam_msg)
        f.write('\n')
    
    #subprocess.run(cmd_classify_beam, shell=True, check=True)
    #subprocess.run(cmd_sim_rad)
