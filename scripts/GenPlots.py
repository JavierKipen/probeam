# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 11:30:36 2023

@author: JK-WORK
"""
import numpy as np
import csv
import pandas as pd
from sklearn import metrics
import time
import matplotlib.pyplot as plt

datasets_path="C:/Users/JK-WORK/Documents/modifWhatprot/Own/HMM_modif/Datasets/ForPaper/";

n_proteins=[20,50,100,200,500,1000,2000,5000,10000,20000];

def get_crossval_acc(true_ids, y_pred,n_folds=10):
    folds_ids=np.array_split(true_ids, n_folds)
    folds_y_pred=np.array_split(y_pred, n_folds)
    accs=[]
    for i in range(n_folds):
        accs.append(metrics.accuracy_score(folds_ids[i], folds_y_pred[i]))
    return np.mean(accs),np.std(accs)/np.sqrt(n_folds)

def get_timing_res(path,n_samples):
    f = open(path, "r")
    msg_hybrid=f.readline();
    msg_beam=f.readline();
    time_hybrid = float(msg_hybrid[msg_hybrid.find(":")+2:-1])/n_samples
    time_beam = float(msg_beam[msg_beam.find(":")+2:-1])*1e-6
    return time_beam,time_hybrid

mean_accs_beam=[]; mean_accs_beam_std=[];
mean_accs_hybrid=[]; mean_accs_hybrid_std=[];
times_beam=[];times_hybrid=[];
for n_prot in n_proteins:
    curr_ds=datasets_path+str(n_prot)+"Prot/"
    true_ids =  pd.read_csv(curr_ds+'true-ids.tsv', sep='\t').to_numpy().flatten()
    predictions_beam = pd.read_csv(curr_ds+"BeamSearchPred10.csv")[' best_pep_iz'].to_numpy();
    predictions_hybrid = pd.read_csv(curr_ds+"predictionsHybrid.csv")['best_pep_iz'].to_numpy();
    mean_acc_beam, mean_acc_beam_std = get_crossval_acc(true_ids,predictions_beam)
    mean_acc_hybrid, mean_acc_hybrid_std = get_crossval_acc(true_ids,predictions_hybrid)
    mean_accs_beam.append(mean_acc_beam);mean_accs_beam_std.append(mean_acc_beam_std)
    mean_accs_hybrid.append(mean_acc_hybrid);mean_accs_hybrid_std.append(mean_acc_hybrid_std)
    time_beam,time_hybrid=get_timing_res(curr_ds+"timing-results.txt",len(true_ids))
    times_beam.append(time_beam);times_hybrid.append(time_hybrid)
    
    
plt.figure(dpi=300)
plt.errorbar(n_proteins,mean_accs_beam,yerr=mean_accs_beam_std,label="Beam")
plt.errorbar(n_proteins,mean_accs_hybrid,yerr=mean_accs_hybrid_std,label="Whatprot")
plt.legend()
plt.xscale("log")
plt.xlabel("Number of proteins")
plt.ylabel("Accuracy")
plt.grid()
plt.savefig("../Accuracy.png")
plt.show()

plt.figure(dpi=300)
plt.loglog(n_proteins,times_beam,label="Beam")
plt.loglog(n_proteins,times_hybrid,label="Whatprot")
plt.legend()
plt.xlabel("Number of proteins")
plt.ylabel("Run time for a single read [s]")
plt.grid()
plt.savefig("../Runtime.png")
plt.show()