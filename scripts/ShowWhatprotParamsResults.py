# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 13:35:44 2023

@author: JK-WORK
"""

import numpy as np
import csv
import pandas as pd
from sklearn import metrics
import time
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os 
import ipdb

#datasets_path="C:/Users/JK-WORK/Documents/modifWhatprot/Own/HMM_modif/Datasets/ForPaper/";
#datasets_path="C:/Users/JK-WORK/Desktop/probeam/probeam/data/NormDatasets/"
#datasets_path="../data/NormDatasets/"
datasets_path="../data/LongDatasets/"
n_proteins=20000;

def get_crossval_acc(true_ids, y_pred, dye_seqs_map, n_folds=10):
    folds_ids=np.array_split(true_ids, n_folds)
    folds_y_pred=np.array_split(y_pred, n_folds)
    accs=[]
    for i in range(n_folds):
        total_acc = 0
        for j in range(len(folds_ids[i])):
            if folds_ids[i][j] == folds_y_pred[i][j]:
                total_acc = total_acc + (1/dye_seqs_map[folds_ids[i][j]])
        total_acc = total_acc/len(folds_ids[i]);
        accs.append(total_acc)
    return np.mean(accs),np.std(accs)/np.sqrt(n_folds)

def get_param(file_path,param,sep_ini="_",sep_end="_"):
    sub_str=file_path[file_path.find(param+sep_ini)+2:];
    if sub_str.find(sep_end)==-1: #If doesnt find final character, assumes a dot too.
        param_value_str=sub_str[:sub_str.find(".")]
    else:
        param_value_str=sub_str[:sub_str.find(sep_end)]
    param=int(param_value_str);
    return param

def get_map_times(timing_path,n_samples):
    f = open(timing_path, "r")
    map_out={};
    cont_read=True;
    while cont_read: 
        line=f.readline();
        if line:
            sub_strk=line[line.find("K=")+2:];
            k_str=sub_strk[:sub_strk.find(" ")]
            sub_strh=line[line.find("H= ")+3:];
            h_str=sub_strh[:sub_strh.find(":")]
            t_str=sub_strh[sub_strh.find(":")+2:sub_strh.find('\n')]
            t=float(t_str)/n_samples;k=int(k_str);h=int(h_str)
            map_out[tuple([k,h])]=t;
        else:
            cont_read=False;
    return map_out

curr_ds=datasets_path+str(n_proteins)+"Prot/"

dye_seqs=pd.read_csv(curr_ds+"dye-seqs.tsv", sep='\t', skiprows=2,header=None)

dye_seqs_map = {dye_seqs.loc[i,2]:dye_seqs.loc[i,1] for i in range(len(dye_seqs))} #Maps pep id to count of peptides same pep id

count=0;
#ipdb.set_trace();
ks=[];hs=[];Accs=[];std_accs=[]; ts=[];
#ipdb.set_trace();
true_ids =  pd.read_csv(curr_ds+'true-ids.tsv', sep='\t').to_numpy().flatten()
n_samples=len(true_ids);
timing_map=get_map_times(curr_ds+'pWtiming-results.txt',n_samples);

print("Showing results from  folder " + curr_ds);
for file_path in os.listdir(curr_ds):
    if file_path.startswith("predictionsHybrid_"):
        k=get_param(file_path,'K');
        h=get_param(file_path,'H');
        predictions_hybrid = pd.read_csv(curr_ds+file_path)['best_pep_iz'].to_numpy();
        acc_hybrid,std_hybrid=get_crossval_acc(true_ids, predictions_hybrid, dye_seqs_map);
        if tuple([k,h]) in timing_map:
            ts.append(timing_map[tuple([k,h])]); #Now not needed
        else:
            ts.append(0); #Now not needed
        ks.append(k);hs.append(h);Accs.append(acc_hybrid);std_accs.append(std_hybrid);

dict= {'K': ks, 'H' : hs, 'Acc': Accs,'Acc std': std_accs, 't': ts}
df = pd.DataFrame(dict)

print(df.sort_values(by=["t","K","H"]))
