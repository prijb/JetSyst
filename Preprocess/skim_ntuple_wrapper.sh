#!/bin/bash
#This wrapper runs the python script for running on a file
#Go to the CMSSW directory and set up the environment
cd $1
source /cvmfs/cms.cern.ch/cmsset_default.sh
source /vols/grid/cms/setup.sh
export X509_USER_PROXY=$1/condor_submission/cms.proxy
cmsenv
#Run the python script
python3 Preprocess/skim_ntuple_muonjet.py -i $2 -o $3
