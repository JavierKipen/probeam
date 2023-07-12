# ProBeam 
A beam-search decoder for peptide classification in fluorosequencing by Javier Kipen and Prof. Joakim Jaldén

Kungliga Tekniska Högskolan (KTH), Stockholm, Sweden

## Reproducing results

The steps to reproduce the results obtained in this paper are listed below. Anaconda is a requirement to run the scripts ([Anaconda Installation](https://docs.anaconda.com/free/anaconda/install/linux/)), and also LaTeX must be installed for the plots. The last one can be installed with the following command:
```
sudo apt-get install dvipng texlive-latex-extra texlive-fonts-recommended cm-super
```

First, open a terminal and enter the probeam folder. An anaconda environment should be created and activated with the following code:

```
conda env create -f env/environment.yml --n probeam
conda activate probeam
```
Then the probeam code has to be compiled, which is done by

```
cd code/
make probeam
cd ..
```
The following script clones whatprot and builds it:
```
bash scripts/setupWhatprot
```

Finally the codes to create datasets, classify the reads and plot the results are python scripts. All these scripts should be run from the probeam folder:
```
python3 scripts/DatasetGen.py %Generates all datasets for the results
python3 scripts/DatasetClassify.py %Runs probeam and whatprot on the datasets from different numbers of proteins
python3 scripts/Genplots.py %Processes the classifications obtained and generates the plots of the papers
python3 scripts/ForHMMComparison.py %Runs the classifications needed to compare with the HMM (without kNN prefilter)

python3 scripts/BestNbSearch.py %Different number of beans run on the 20k Prot dataset
python3 scripts/analyzeNB.py %Generates a csv from the results
python3 scripts/GenPlots.py %Generates the figures of the paper with the results
```
# Reusing the code

Upon use of this code for research purposes, please cite [our paper](https://docs.anaconda.com/free/anaconda/install/linux/)
