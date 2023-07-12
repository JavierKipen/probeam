import numpy as np
import csv
import pandas as pd
from sklearn import metrics
import time
#import matplotlib.pyplot as plt

path_datasets="data/NormDatasets/"
n_proteins=20000;
path_big_ds= path_datasets + str(n_proteins)+"Prot/"

def get_crossval_acc(true_ids, y_pred,n_folds=10):
    folds_ids=np.array_split(true_ids, n_folds)
    folds_y_pred=np.array_split(y_pred, n_folds)
    accs=[]
    for i in range(n_folds):
        accs.append(metrics.accuracy_score(folds_ids[i], folds_y_pred[i]))
    return np.mean(accs),np.std(accs)/np.sqrt(n_folds)


df=pd.read_csv(path_big_ds+"NBtiming-results.txt", header=None)
df['Accuracy'] = 0;
df.columns=["N beam", "Microsecs per read", "Accuracy"];
for i in range(len(df)):
    true_ids =  pd.read_csv(path_big_ds+'true-ids.tsv', sep='\t').to_numpy().flatten()
    predictions_beam = pd.read_csv(path_big_ds+"BeamSearchPred"+str(int(df.iloc[i]["N beam"]))+".csv")[' best_pep_iz'].to_numpy();
    mean_acc_beam, mean_acc_beam_std = get_crossval_acc(true_ids,predictions_beam)
    df.loc[i,"Accuracy"]=mean_acc_beam;

print(df)

df.to_csv("results/NbTable20000Prot.csv",index=False)
