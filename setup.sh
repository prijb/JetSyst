#!/bin/bash
CONDOR_SUBMISSION_DIR="condor_submission"
OUTPUT_DIR_NAME="outputs"
PLOT_DIR_NAME="plots"
#Setup output directory
#Make them in the current working directory if they don't exist
if [ ! -d "${OUTPUT_DIR_NAME}" ]; then
    mkdir ${OUTPUT_DIR_NAME}
fi
if [ ! -d "${PLOT_DIR_NAME}" ]; then
    mkdir ${PLOT_DIR_NAME}
fi
if [ ! -d "${CONDOR_SUBMISSION_DIR}" ]; then
    mkdir ${CONDOR_SUBMISSION_DIR}
fi
#CMSSW setup
if [ -z ${CMSSW_BASE+x} ]; then
    cmsenv
    echo "CMSSW environment set up in ${CMSSW_BASE}"
else
    echo "CMSSW environment already set up in ${CMSSW_BASE}"
fi
#Check if proxy is set up
if [ -z ${X509_USER_PROXY+x} ]; then
    echo "Setting up proxy"
#    voms-proxy-init --rfc --voms cms --valid 192:00
    voms-proxy-init --rfc --voms cms --valid 192:00 --out ${HOME}/cms.proxy; export X509_USER_PROXY=${HOME}/cms.proxy
else
    echo "Proxy already set up"
fi
#Copy proxy to local condor submission directory
cp ${X509_USER_PROXY} condor_submission/cms.proxy
