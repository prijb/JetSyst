# Run fit jobs for all pT and eta bins
import os
import argparse
import numpy as np

# Encoding function to convert -3.0 to "neg3p0" and 3.0 to "pos3p0"
def encode_eta(eta):
    if eta < 0:
        return "neg" + str(abs(eta)).replace(".", "p")
    else:
        return "pos" + str(eta).replace(".", "p")

# Argument parser
# Example: python3 Preprocess/fit_histos.py --input outputs/Histogram/muon2024g_muonjet --output outputs/Fit/muon2024g_muonjet --script Preprocess/fit_histo.py --plotdir plots/Fit/muon2024g_muonjet
parser = argparse.ArgumentParser(description='Create histograms')
parser.add_argument('--input', '-i', type=str, help='Input directory')
parser.add_argument('--output', '-o', type=str, help='Output directory')
parser.add_argument('--plotdir', '-p', type=str, default="plots/Fit/data", help='Plot directory')
parser.add_argument('--runlocal', action='store_true', help='Run locally without condor')
parser.add_argument("--script", '-sc', type=str, default="Preprocess/fit_histo.py", help="Histogramming script")

args = parser.parse_args()
cwd = os.getcwd()
script = args.script
pt_edges = np.concatenate([np.arange(20, 60, 4), np.arange(60, 100, 5), np.arange(100, 200, 10), np.arange(200, 300, 20), np.arange(300, 500, 50)])
eta_edges = np.array([-3.0, -2.5, -2.0, -1.3, -0.5, 0.0, 0.5, 1.3, 2.0, 2.5, 3.0])
num_pt_bins = len(pt_edges) - 1
num_eta_bins = len(eta_edges) - 1

# Make a directory for the output path
os.makedirs(args.output, exist_ok=True)
os.makedirs(args.plotdir, exist_ok=True)

#Directories and files for condor submission
condor_parent_dir = "condor_submission/Fit"
os.makedirs(f"{condor_parent_dir}", exist_ok=True)
os.makedirs(f"{condor_parent_dir}/logs", exist_ok=True)

# Make the input arguments file for the condor job
with open(f"{condor_parent_dir}/fit_histos_args.txt", "w") as f:
    for i_eta_bin in range(len(eta_edges)-1):
        eta_min = eta_edges[i_eta_bin]
        eta_max = eta_edges[i_eta_bin+1]

        for i_pt_bin in range(len(pt_edges)-1):
            pt_min = pt_edges[i_pt_bin]
            pt_max = pt_edges[i_pt_bin+1]

            outfile_name = f"fit_pt_{i_pt_bin}_eta_{i_eta_bin}.pkl"
            f.write(f"{args.input}/histos_pt_{i_pt_bin}_eta_{i_eta_bin}.root {args.output}/{outfile_name} {int(pt_min)} {int(pt_max)} {encode_eta(eta_min)} {encode_eta(eta_max)} 1 {args.plotdir}\n")

        # pT inclusive jobs
        outfile_name = f"fit_pt_inclusive_eta_{i_eta_bin}.pkl"
        f.write(f"{args.input}/histos_pt_inclusive_eta_{i_eta_bin}.root {args.output}/{outfile_name} {int(pt_edges[0])} {int(pt_edges[-1])} {encode_eta(eta_min)} {encode_eta(eta_max)} 1 {args.plotdir}\n")

    # eta inclusive jobs
    for i_pt_bin in range(len(pt_edges)-1):
        pt_min = pt_edges[i_pt_bin]
        pt_max = pt_edges[i_pt_bin+1]

        outfile_name = f"fit_eta_inclusive_pt_{i_pt_bin}.pkl"
        f.write(f"{args.input}/histos_eta_inclusive_pt_{i_pt_bin}.root {args.output}/{outfile_name} {int(pt_min)} {int(pt_max)} {encode_eta(eta_edges[0])} {encode_eta(eta_edges[-1])} 1 {args.plotdir}\n")

# Create the wrapper file
wrapper_file_content = f"""#!/bin/bash
#This wrapper runs the python script for running on a file
#Go to the CMSSW directory and set up the environment
cd {cwd}
source /cvmfs/cms.cern.ch/cmsset_default.sh
source /vols/grid/cms/setup.sh
export X509_USER_PROXY={cwd}/condor_submission/cms.proxy
cmsenv
#Run the python script
python3 {script} --input $1 --output $2 --ptlow $3 --pthigh $4 --etalow $5 --etahigh $6 --rebin $7 --plotdir $8
"""

with open("Preprocess/fit_histo_wrapper.sh", "w") as f:
    f.write(wrapper_file_content)
os.system("chmod +x Preprocess/fit_histo_wrapper.sh")

# Create the HTCondor submit file 
submit_file_content = f"""universe = vanilla
executable = {cwd}/Preprocess/fit_histo_wrapper.sh
arguments = $(infile) $(outfile) $(ptlow) $(pthigh) $(etalow) $(etahigh) $(rebin) $(plotdir)
output = {condor_parent_dir}/logs/fit_histo_$(Cluster)_$(Process).out
error = {condor_parent_dir}/logs/fit_histo_$(Cluster)_$(Process).err
log = {condor_parent_dir}/logs/fit_histo_$(Cluster)_$(Process).log
request_cpus = 1
request_memory = 4GB
use_x509userproxy = true
+MaxRuntime = 7199
queue infile, outfile, ptlow, pthigh, etalow, etahigh, rebin, plotdir from {condor_parent_dir}/fit_histos_args.txt
"""

with open(f"{condor_parent_dir}/fit_histos.submit", "w") as f:
    f.write(submit_file_content)

# Delete existing log files
os.system(f"rm {condor_parent_dir}/logs/*")

if args.runlocal:
    print(f"Running jobs locally...")
    with open(f"{condor_parent_dir}/fit_histos_args.txt", "r") as f:
        lines = f.readlines()
        for line in lines:
            parts = line.strip().split()
            infile = parts[0]
            outfile = parts[1]
            ptlow = parts[2]
            pthigh = parts[3]
            etalow = parts[4]
            etahigh = parts[5]
            rebin = parts[6]
            plotdir = parts[7]
            cmd = f"python3 {script} --input {infile} --output {outfile} --ptlow {ptlow} --pthigh {pthigh} --etalow {etalow} --etahigh {etahigh} --rebin {rebin} --plotdir {plotdir}"
            print(f"Running: {cmd}")
            os.system(f"{cmd}")

# Run the condor job
else:
    os.system(f"condor_submit {condor_parent_dir}/fit_histos.submit")