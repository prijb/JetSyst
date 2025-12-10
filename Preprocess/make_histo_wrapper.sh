#!/bin/bash
#This wrapper runs the python script for running on a file
#Go to the CMSSW directory and set up the environment
cd /vols/cms/pb4918/L1Scouting/Dec25/JetSyst/CMSSW_14_1_0_pre4/src/JetSyst
source /cvmfs/cms.cern.ch/cmsset_default.sh
source /vols/grid/cms/setup.sh
export X509_USER_PROXY=$1/condor_submission/cms.proxy
cmsenv
#Run the python script
python3 Preprocess/make_histo.py --input $1 --output $2 --ptrange $3 $4 --etarange $5 $6
