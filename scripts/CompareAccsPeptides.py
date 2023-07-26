# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 13:35:44 2023

@author: JK-WORK
"""

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
from matplotlib.colors import LogNorm

#datasets_path="C:/Users/JK-WORK/Documents/modifWhatprot/Own/HMM_modif/Datasets/ForPaper/";
#datasets_path="C:/Users/JK-WORK/Desktop/probeam/probeam/data/NormDatasets/"
datasets_path="../data/NormDatasets/"
n_proteins=1000;

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


curr_ds=datasets_path+str(n_proteins)+"Prot/"

dye_seqs=pd.read_csv(curr_ds+"dye-seqs.tsv", sep='\t', skiprows=2,header=None)

dye_seqs_map = {dye_seqs.loc[i,2]:dye_seqs.loc[i,1] for i in range(len(dye_seqs))} #Maps pep id to count of peptides same pep id

predictions_beam = pd.read_csv(curr_ds+"BeamSearchPred15.csv")[' best_pep_iz'].to_numpy();
#predictions_beam_old = pd.read_csv(curr_ds+"BeamSearchPred15Old.csv")[' best_pep_iz'].to_numpy();
predictions_hybrid = pd.read_csv(curr_ds+"predictionsHybrid.csv")['best_pep_iz'].to_numpy();
predictions_hmm = pd.read_csv(curr_ds+"predictionsHMM.csv")['best_pep_iz'].to_numpy();
#predictions_kNN = pd.read_csv(curr_ds+"predictionskNN.csv")['best_pep_iz'].to_numpy();

true_ids =  pd.read_csv(curr_ds+'true-ids.tsv', sep='\t').to_numpy().flatten()

acc_beam,std_beam=get_crossval_acc(true_ids, predictions_beam, dye_seqs_map);
#acc_beam_old,std_beam_old=get_crossval_acc(true_ids, predictions_beam_old, dye_seqs_map);
acc_hybrid,std_hybrid=get_crossval_acc(true_ids, predictions_hybrid, dye_seqs_map);
acc_hmm,std_hmm=get_crossval_acc(true_ids, predictions_hmm, dye_seqs_map);
#acc_knn,std_knn=get_crossval_acc(true_ids, predictions_kNN, dye_seqs_map);

print(str(acc_beam) + str(std_beam))
#print(str(acc_beam_old) + str(std_beam_old))
print(str(acc_hybrid) + str(std_hybrid))
print(str(acc_hmm) + str(std_hmm))
#print(str(acc_knn) + str(std_knn))
