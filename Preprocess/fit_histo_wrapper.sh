#!/bin/bash
#This wrapper runs the python script for running on a file
#Go to the CMSSW directory and set up the environment
cd /vols/cms/pb4918/L1Scouting/Dec25/JetSyst/CMSSW_14_1_0_pre4/src/JetSyst
source /cvmfs/cms.cern.ch/cmsset_default.sh
source /vols/grid/cms/setup.sh
export X509_USER_PROXY=/vols/cms/pb4918/L1Scouting/Dec25/JetSyst/CMSSW_14_1_0_pre4/src/JetSyst/condor_submission/cms.proxy
cmsenv
#Run the python script
python3 Preprocess/fit_histo.py --input $1 --output $2 --ptlow $3 --pthigh $4 --etalow $5 --etahigh $6 --rebin $7 --plotdir $8
