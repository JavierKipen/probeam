import os, shutil
import subprocess

##Script to have the data for figures 3 and 4 (Comparison of score HMM vs beam dec, and comparison of HMM accuracy vs Beam accuracy)

path_data=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+"/data/";
path_data_common=path_data +"common/"
path_data_norm=path_data+"NormDatasets/";
path_data_long=path_data+"LongDatasets/";
path_seq_params = path_data_common + "seq-params.json"

#HMM run for 1000Prot to compare scores
path_score_ds=path_data_norm+"1000Prot/";
path_dye_seqs= path_score_ds +"dye-seqs.tsv";
path_radiometries= path_score_ds +"radiometries.tsv";
path_output=path_score_ds+"predictionsHMM.csv"
cmd_classify_HMM = "./bin/whatprot classify hmm -P " + path_seq_params + " -S " + path_dye_seqs + " -R " + path_radiometries + " -Y " + path_output;
output=subprocess.check_output(cmd_classify_HMM, shell=True)
print(str(output))



#To compare accuracies!
path_score_ds=path_data_long+"1000Prot/";
path_dye_seqs= path_score_ds +"dye-seqs.tsv";
path_radiometries= path_score_ds +"radiometries.tsv";
path_output=path_score_ds+"predictionsHMM.csv"
cmd_classify_HMM = "./bin/whatprot classify hmm -P " + path_seq_params + " -S " + path_dye_seqs + " -R " + path_radiometries + " -Y " + path_output;
output=subprocess.check_output(cmd_classify_HMM, shell=True)
print(str(output))
cmd_classify_beam = "./bin/probeam " + path_score_ds + " -b 15"; #NB=15
output=subprocess.check_output(cmd_classify_beam, shell=True)
print(str(output))
