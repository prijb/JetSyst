# Run histogram creation jobs for all pT and eta bins
import os
import argparse
import numpy as np

#Argument parser
# Example: python3 Preprocess/make_histos.py --input outputs/Skim/merged/muon2024g_muonjet.root --output outputs/Histogram/muon2024g_muonjet --script Preprocess/make_histo.py
parser = argparse.ArgumentParser(description='Create histograms')
parser.add_argument('--input', '-i', type=str, help='Input file')
parser.add_argument('--output', '-o', type=str, help='Output directory')
parser.add_argument("--script", '-sc', type=str, default="Preprocess/make_histo.py", help="Histogramming script")

args = parser.parse_args()
cwd = os.getcwd()
script = args.script
pt_edges = np.concatenate([np.arange(20, 60, 4), np.arange(60, 100, 5), np.arange(100, 200, 10), np.arange(200, 300, 20), np.arange(300, 500, 50)])
eta_edges = np.array([-3.0, -2.5, -2.0, -1.3, -0.5, 0.0, 0.5, 1.3, 2.0, 2.5, 3.0])

#Make a directory for the output path
os.makedirs(args.output, exist_ok=True)

#Directories and files for condor submission
condor_parent_dir = "condor_submission/Histogram"
os.makedirs(f"{condor_parent_dir}", exist_ok=True)
os.makedirs(f"{condor_parent_dir}/logs", exist_ok=True)

#Make the input arguments file for the condor job
with open(f"{condor_parent_dir}/make_histos_args.txt", "w") as f:
    for i_eta_bin in range(len(eta_edges)-1):
        eta_min = eta_edges[i_eta_bin]
        eta_max = eta_edges[i_eta_bin+1]

        for i_pt_bin in range(len(pt_edges)-1):
            pt_min = pt_edges[i_pt_bin]
            pt_max = pt_edges[i_pt_bin+1]

            outfile_name = f"histos_pt_{i_pt_bin}_eta_{i_eta_bin}.root"
            f.write(f"{args.input} {args.output}/{outfile_name} {pt_min} {pt_max} {eta_min} {eta_max}\n")
    
    # pT inclusive jobs
    for i_eta_bin in range(len(eta_edges)-1):
        eta_min = eta_edges[i_eta_bin]
        eta_max = eta_edges[i_eta_bin+1]

        outfile_name = f"histos_pt_inclusive_eta_{i_eta_bin}.root"
        f.write(f"{args.input} {args.output}/{outfile_name} 0 200 {eta_min} {eta_max}\n")

    # eta inclusive jobs
    for i_pt_bin in range(len(pt_edges)-1):
        pt_min = pt_edges[i_pt_bin]
        pt_max = pt_edges[i_pt_bin+1]

        outfile_name = f"histos_eta_inclusive_pt_{i_pt_bin}.root"
        f.write(f"{args.input} {args.output}/{outfile_name} {pt_min} {pt_max} -3.0 3.0\n")

#Create the wrapper file
wrapper_file_content = f"""#!/bin/bash
#This wrapper runs the python script for running on a file
#Go to the CMSSW directory and set up the environment
cd {cwd}
source /cvmfs/cms.cern.ch/cmsset_default.sh
source /vols/grid/cms/setup.sh
export X509_USER_PROXY=$1/condor_submission/cms.proxy
cmsenv
#Run the python script
python3 {script} --input $1 --output $2 --ptrange $3 $4 --etarange $5 $6
"""

with open("Preprocess/make_histo_wrapper.sh", "w") as f:
    f.write(wrapper_file_content)
os.system("chmod +x Preprocess/make_histo_wrapper.sh")

#Create the HTCondor submit file
submit_file_content = f"""universe = vanilla
executable = {cwd}/Preprocess/make_histo_wrapper.sh
arguments = $(infile) $(outfile) $(pt_min) $(pt_max) $(eta_min) $(eta_max)
output = {condor_parent_dir}/logs/make_histo_$(Cluster)_$(Process).out
error = {condor_parent_dir}/logs/make_histo_$(Cluster)_$(Process).err
log = {condor_parent_dir}/logs/make_histo_$(Cluster)_$(Process).log
request_cpus = 1
request_memory = 4GB
use_x509userproxy = true
+MaxRuntime = 7199
queue infile, outfile, pt_min, pt_max, eta_min, eta_max from {condor_parent_dir}/make_histos_args.txt
"""

with open(f"{condor_parent_dir}/make_histos.submit", "w") as f:
    f.write(submit_file_content)

#Delete existing log files
os.system(f"rm {condor_parent_dir}/logs/*")

#Run the condor job
os.system(f"condor_submit {condor_parent_dir}/make_histos.submit")