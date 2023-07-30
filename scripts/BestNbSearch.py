# -*- coding: utf-8 -*-
"""
Created on Fri May 26 15:19:30 2023

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

n_proteins=20000; #Only on the 20k dataset
n_beams=[2,3,5,7,10,12,15,20,40];

path_datasets_figure=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+"/data/NormDatasets/";

path_datasets_figure_common_files= path_datasets_figure +"CommonFiles/"
prot_fasta=path_datasets_figure_common_files+"UP000005640_9606.fasta"
seq_params_path= path_datasets_figure_common_files + "seq-params.json"


protein_folder=path_datasets_figure+str(n_proteins)+"Prot/"

whole_message=""
for n_beam in n_beams:
    print("Number of beams: "+ str(n_beams))
    
    dye_seqs_path=protein_folder+"dye-seqs.tsv";
    dye_tracks_path=protein_folder+"dye-tracks.tsv";
    radiometries_path=protein_folder+"radiometries.tsv";
    true_ids_path=protein_folder+"true-ids.tsv";
    predictions_hybrid_path = protein_folder+"predictionsHybrid.csv"
    
    cmd_classify_beam = "./bin/probeam " + path_datasets_figure + str(n_proteins)+"Prot/"+ " -b "+str(n_beam);
    
    #output=subprocess.run(cmd_classify_hybrid, shell=True, check=True)
    outputBeam=subprocess.check_output(cmd_classify_beam, shell=True)
    #outputBeam="Time: 201 microseconds"
    computing_time_beam=get_proc_time_beam(outputBeam)
    beam_msg=str(n_beam) + ","+ str(computing_time_beam)+ "\n";
    print(str(outputBeam))
    whole_message = whole_message+beam_msg;
    
with open(protein_folder+"NBtiming-results.txt", 'w') as f:
    f.write(whole_message)
    f.write('\n')

    #subprocess.run(cmd_classify_beam, shell=True, check=True)
    #subprocess.run(cmd_sim_rad)
