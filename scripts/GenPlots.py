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
import ipdb
from matplotlib.colors import LogNorm

#datasets_path="C:/Users/JK-WORK/Documents/modifWhatprot/Own/HMM_modif/Datasets/ForPaper/";
datasets_path="data/NormDatasets/"

n_proteins=[20,50,100,200,500,1000,2000,5000,10000,20000];

def get_crossval_acc(true_ids, y_pred,dye_seqs_map,n_folds=10):
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
        #accs.append(metrics.accuracy_score(folds_ids[i], folds_y_pred[i]))
    return np.mean(accs),np.std(accs)/np.sqrt(n_folds)

###Figure 1 ####

def get_timing_res(path,n_samples):
    f = open(path, "r")
    msg_hybrid=f.readline();
    msg_kNN=f.readline();
    msg_beam_1=f.readline();
    msg_beam_2=f.readline();
    nb1=int(msg_beam_1[:1])
    nb2=int(msg_beam_2[:2])
    time_hybrid = float(msg_hybrid[msg_hybrid.find(":")+2:-1])/n_samples
    time_knn = float(msg_kNN[msg_kNN.find(":")+2:-1])/n_samples
    time_beam_1 = float(msg_beam_1[msg_beam_1.find(":")+2:-1])*1e-6
    time_beam_2 = float(msg_beam_2[msg_beam_2.find(":")+2:-1])*1e-6
    return time_beam_1,time_beam_2,time_hybrid,time_knn,nb1,nb2

mean_accs_beam_1=[]; mean_accs_beam_std_1=[];
mean_accs_beam_2=[]; mean_accs_beam_std_2=[];
mean_accs_hybrid=[]; mean_accs_hybrid_std=[];
mean_accs_knn=[]; mean_accs_knn_std=[];
times_beam_1=[];times_beam_2=[];times_hybrid=[];times_knn=[];
for n_prot in n_proteins:
    curr_ds=datasets_path+str(n_prot)+"Prot/"
    true_ids =  pd.read_csv(curr_ds+'true-ids.tsv', sep='\t').to_numpy().flatten()
    time_beam_1,time_beam_2,time_hybrid,time_knn,nb1,nb2=get_timing_res(curr_ds+"timing-results.txt",len(true_ids));
    predictions_beam_1 = pd.read_csv(curr_ds+"BeamSearchPred"+str(nb1)+".csv")[' best_pep_iz'].to_numpy();
    predictions_beam_2 = pd.read_csv(curr_ds+"BeamSearchPred"+str(nb2)+".csv")[' best_pep_iz'].to_numpy();
    predictions_hybrid = pd.read_csv(curr_ds+"predictionsHybrid.csv")['best_pep_iz'].to_numpy();
    predictions_knn = pd.read_csv(curr_ds+"predictionskNN.csv")['best_pep_iz'].to_numpy();
    
    dye_seqs=pd.read_csv(curr_ds+"dye-seqs.tsv", sep='\t', skiprows=2,header=None)
    dye_seqs_map = {dye_seqs.loc[i,2]:dye_seqs.loc[i,1] for i in range(len(dye_seqs))} #Maps pep id to count of peptides same pep id

    mean_acc_beam_1, mean_acc_beam_std_1 = get_crossval_acc(true_ids,predictions_beam_1,dye_seqs_map)
    mean_acc_beam_2, mean_acc_beam_std_2 = get_crossval_acc(true_ids,predictions_beam_2,dye_seqs_map)
    mean_acc_hybrid, mean_acc_hybrid_std = get_crossval_acc(true_ids,predictions_hybrid,dye_seqs_map)
    mean_acc_knn, mean_acc_knn_std = get_crossval_acc(true_ids,predictions_knn,dye_seqs_map)
    mean_accs_beam_1.append(mean_acc_beam_1);mean_accs_beam_std_1.append(mean_acc_beam_std_1)
    mean_accs_beam_2.append(mean_acc_beam_2);mean_accs_beam_std_2.append(mean_acc_beam_std_2)
    mean_accs_hybrid.append(mean_acc_hybrid);mean_accs_hybrid_std.append(mean_acc_hybrid_std)
    mean_accs_knn.append(mean_acc_knn);mean_accs_knn_std.append(mean_acc_knn_std)
    times_beam_1.append(time_beam_1);times_beam_2.append(time_beam_2);times_hybrid.append(time_hybrid);times_knn.append(time_knn);
    
plt.rcParams.update({
    "text.usetex": True})
plt.figure(dpi=300)
plt.errorbar(n_proteins,mean_accs_beam_1,yerr=mean_accs_beam_std_1,label=r'Beam $\textrm{N}_{\textrm{B}}=7$')
plt.errorbar(n_proteins,mean_accs_beam_2,yerr=mean_accs_beam_std_2,label=r'Beam $\textrm{N}_{\textrm{B}}=20$')
plt.errorbar(n_proteins,mean_accs_hybrid,yerr=mean_accs_hybrid_std,label="Whatprot")
plt.errorbar(n_proteins,mean_accs_knn,yerr=mean_accs_knn_std,label="kNN")
plt.legend()
plt.xscale("log")
plt.xlabel("Number of proteins")
plt.ylabel("Accuracy")
plt.grid()
plt.savefig("results/Accuracy.png")
#plt.show()

plt.figure(dpi=300)
plt.loglog(n_proteins,times_beam_1,label=r'Beam $\textrm{N}_{\textrm{B}}=7$')
plt.loglog(n_proteins,times_beam_2,label=r'Beam $\textrm{N}_{\textrm{B}}=20$')
plt.loglog(n_proteins,times_hybrid,label="Whatprot")
plt.loglog(n_proteins,times_knn,label="kNN")
plt.legend()
plt.xlabel("Number of proteins")
plt.ylabel("Run time for a single read [s]")
plt.grid()
plt.savefig("results/Runtime.png")

###Figure 2 ####
#ipdb.set_trace()
"""
dataset_thousand_path=datasets_path+"1000Prot/";

true_ids =  pd.read_csv(dataset_thousand_path+'true-ids.tsv', sep='\t').to_numpy().flatten()
HMM = pd.read_csv(dataset_thousand_path+'predictionsHMM.csv')
y_predHMM=HMM['best_pep_iz'].to_numpy(); 
y_predHMMScore=HMM['best_pep_score'].to_numpy(); 
Beam = pd.read_csv(dataset_thousand_path+'BeamSearchPred15.csv')
y_predB=Beam[' best_pep_iz'].to_numpy(); 
y_predBScore=Beam[' best_pep_score'].to_numpy(); 

def doHist2D(ypredX,ypredY,varX,varY,labelX,labelY,ids):
    len_2_use=min(len(varX),len(varY));
    plt.figure(dpi=300)
    x_vals=[];y_vals=[];
    for i in range(len_2_use):
        if ypredX[i]==ypredY[i] and ypredX[i]==ids[i]:
            x_vals.append(varX[i])
            y_vals.append(varY[i])
    h=plt.hist2d(x_vals,y_vals, norm=LogNorm(),bins=50)
    plt.colorbar(h[3])
    plt.xlabel(labelX)
    plt.ylabel(labelY)

doHist2D(y_predB,y_predHMM,y_predBScore,y_predHMMScore,r'Beam $\textrm{N}_{\textrm{B}}=15$',"HMM",true_ids)
plt.savefig("results/histComp.eps")
"""
##Accuracy for table 

"""

dataset_long="data/LongDatasets/1000Prot/"
#dataset_long="data/NormDatasets/1000Prot/"
predictions_beam = pd.read_csv(dataset_long+"BeamSearchPred15.csv")[' best_pep_iz'].to_numpy();
predictions_hmm = pd.read_csv(dataset_long+"predictionsHMM.csv")['best_pep_iz'].to_numpy();
true_ids =  pd.read_csv(dataset_long+'true-ids.tsv', sep='\t').to_numpy().flatten()

mean_acc_beam, mean_acc_beam_std = get_crossval_acc(true_ids,predictions_beam)
mean_acc_hmm, mean_acc_hmm_std = get_crossval_acc(true_ids,predictions_hmm)

print("Results for long dataset: ")
print("Beam Decoder - Accuracy " +  str(mean_acc_beam) + ", Std: " + str(mean_acc_beam_std))
print("HMM - Accuracy " +  str(mean_acc_hmm) + ", Std: " + str(mean_acc_hmm_std))

"""
