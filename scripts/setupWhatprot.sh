#!/bin/bash
#Bash to download whatprot in the biorxiv version, compile it and copy the executables to the bin folder
#It is supposed to be run from the folder inside probeam.
cd ext
git clone https://github.com/marcottelab/whatprot.git
cd whatprot 
git checkout eab71d5
echo "Whatprot downloaded and checkout to corresponding commit"
cd cc_code
make release
cd ..
cp cc_code/bin/release/whatprot ../../bin/whatprot
