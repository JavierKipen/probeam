# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 17:10:28 2023

@author: JK-WORK
"""

import os, shutil
import subprocess
import ipdb

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

n_proteins=20000;

path_data=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+"/data/";
#full_path_to_ds=path_data+"NormDatasets/";
full_path_to_ds=path_data+"LongDatasets/";
path_datasets_figure_common_files= path_data +"common/"
prot_fasta=path_datasets_figure_common_files+"UP000005640_9606.fasta"
seq_params_path= path_datasets_figure_common_files + "seq-params.json"
protein_folder=full_path_to_ds+str(n_proteins)+"Prot/"
dye_seqs_path=protein_folder+"dye-seqs.tsv";
dye_tracks_path=protein_folder+"dye-tracks.tsv";
radiometries_path=protein_folder+"radiometries.tsv";
true_ids_path=protein_folder+"true-ids.tsv";
cmd_set_threads="export OMP_NUM_THREADS=1 \n";
#cmd_set_threads="export OMP_NUM_THREADS=40 \n";

kH_param_map={100:[25,15],75:[25],200:[15]}
#kH_param_map={200:[20,15]}

msgs_list=[];
for k in kH_param_map.keys():
    h_list=kH_param_map[k];
    for h in h_list:
        print("Running with K=" +str(k)+" H= "+ str(h));
        predictions_hybrid_path = protein_folder+"predictionsHybrid_K_"+str(k)+"_H_"+str(h)+".csv";
        cmd_classify_hybrid = "./bin/whatprot classify hybrid -k "+str(k)+" -s 0.5 -H "+str(h)+" -p 5 -P " + seq_params_path + " -S " + dye_seqs_path + " -T " +dye_tracks_path + " -R " + radiometries_path + " -Y " + predictions_hybrid_path;
        #outputHybrid=subprocess.check_output("export OMP_NUM_THREADS=1\n"+cmd_classify_hybrid, shell=True)
        outputHybrid=subprocess.check_output(cmd_set_threads+cmd_classify_hybrid, shell=True) #To run faster to obtain acc results before
        computing_time_hybrid=get_proc_time(outputHybrid)
        hybrid_msg="Hybrid computing total time (seconds) for K=" +str(k)+" H= "+ str(h)+ ": "+ str(computing_time_hybrid);
        print(str(outputHybrid))
        msgs_list.append(str(hybrid_msg));
with open(protein_folder+"pWtiming-results.txt", 'a') as f:
    for i in range(len(msgs_list)):
        f.write(msgs_list[i])
        f.write('\n')
    #subprocess.run(cmd_classify_beam, shell=True, check=True)
    #subprocess.run(cmd_sim_rad)

