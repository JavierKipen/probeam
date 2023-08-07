import numpy as np
import csv
import pandas as pd
from sklearn import metrics
import time
#import matplotlib.pyplot as plt

#path_datasets="data/NormDatasets/"
path_datasets="data/LongDatasets/"
n_proteins=20000;
path_big_ds= path_datasets + str(n_proteins)+"Prot/"

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
        #accs.append(metrics.accuracy_score(folds_ids[i], folds_y_pred[i]))
    return np.mean(accs),np.std(accs)/np.sqrt(n_folds)

dye_seqs=pd.read_csv(path_big_ds+"dye-seqs.tsv", sep='\t', skiprows=2,header=None)
dye_seqs_map = {dye_seqs.loc[i,2]:dye_seqs.loc[i,1] for i in range(len(dye_seqs))} #Maps pep id to count of peptides same pep id
true_ids =  pd.read_csv(path_big_ds+'true-ids.tsv', sep='\t').to_numpy().flatten()
df=pd.read_csv(path_big_ds+"NBtiming-results.txt", header=None)
df['Accuracy'] = 0;
df['Accuracy std'] = 0;
df.columns=["N beam", "Microsecs per read", "Accuracy", "Accuracy std"];
for i in range(len(df)):
    predictions_beam = pd.read_csv(path_big_ds+"BeamSearchPred"+str(int(df.iloc[i]["N beam"]))+".csv")[' best_pep_iz'].to_numpy();
    mean_acc_beam, mean_acc_beam_std = get_crossval_acc(true_ids,predictions_beam,dye_seqs_map)
    df.loc[i,"Accuracy"]=mean_acc_beam;
    df.loc[i,"Accuracy std"]=mean_acc_beam_std;

print("Results in dataset " + path_big_ds)
print(df)

df.to_csv("results/NbTable20000Prot.csv",index=False)
