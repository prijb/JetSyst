# Extract the scale response from the fits, plot and fit a GPR to the scale
# Plot the jet pT scale
import os 
import pickle
import argparse

import uproot
import hist
import pandas as pd
import numpy as np

# Fitting
import scipy as sp
from scipy.optimize import curve_fit
# Fit using GPR
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C, WhiteKernel, Matern

import matplotlib.pyplot as plt
import mplhep as hep
hep.style.use("CMS")

from matplotlib.ticker import FixedLocator, LogLocator, ScalarFormatter, FuncFormatter
def minor_fmt(val, pos):
    v = int(round(val))
    return str(v) if v in minor_labels else ""

# Argument parser
parser = argparse.ArgumentParser(description="Plot the pT scale from fit results")
parser.add_argument("--input", "-i", type=str, help="Input directory of fits")
parser.add_argument("--output", "-o", type=str, help="Output plot directory")
args = parser.parse_args()

input_dir = args.input
output_dir = args.output
cwd = os.getcwd()

os.makedirs(output_dir, exist_ok=True)

# Make sure these are consistent with production
pt_edges = np.concatenate([np.arange(20, 60, 4), np.arange(60, 100, 5), np.arange(100, 200, 10), np.arange(200, 300, 20), np.arange(300, 500, 50)])
eta_edges = np.array([-3.0, -2.5, -2.0, -1.3, -0.5, 0.0, 0.5, 1.3, 2.0, 2.5, 3.0])
num_pt_bins = len(pt_edges) - 1
num_eta_bins = len(eta_edges) - 1

# Make a pandas dataframe to collect the fit results to plot
col_names = ["pt_bin", "eta_bin", "mean_pt", "mean_scale", "mean_scale_err", "mean_scale_err_bootstrap", "mean_scale_fit", "mean_scale_fit_err", "rms_scale", "std_scale", "reco_pt"]
df = pd.DataFrame(columns=col_names)
df_inclusive = pd.DataFrame(columns=col_names)

# Loop over all the eta bins first, then pT bins and fill dataframe
i_row = 0
for i_eta_bin in range(num_eta_bins):
    for i_pt_bin in range(num_pt_bins):
    #for i_pt_bin in range(2, num_pt_bins):
        fname = f"{input_dir}/fit_pt_{i_pt_bin}_eta_{i_eta_bin}.pkl"

        reco_pt_i = 0.5 * (pt_edges[i_pt_bin] + pt_edges[i_pt_bin+1])

        result = None
        if os.path.exists(fname):
            with open(fname, "rb") as f:
                result = pickle.load(f)
                result = result["scale"]
            
        if result is not None:
            result_i = [int(i_pt_bin), int(i_eta_bin), result["x"], result["mean_guess"], result["mean_guess_err"], result["mean_guess_err_bootstrap"], result["mean_fit"], result["mean_fit_err"], result["std_guess"], result["std_fit"], reco_pt_i]
            df.loc[i_row] = result_i

        else:
            result_i = [int(i_pt_bin), int(i_eta_bin), 0, 0, 0, 0, 0, 0, 0, reco_pt_i]
            #df.loc[i_row] = result_i
            
        i_row += 1

# Load the inclusive eta bin results
for i_pt_bin in range(num_pt_bins):
#for i_pt_bin in range(2, num_pt_bins):
    fname = f"{input_dir}/fit_eta_inclusive_pt_{i_pt_bin}.pkl"

    reco_pt_i = 0.5 * (pt_edges[i_pt_bin] + pt_edges[i_pt_bin+1])

    result = None
    if os.path.exists(fname):
        with open(fname, "rb") as f:
            result = pickle.load(f)
            result = result["scale"]

    if result is not None:
        result_i = [int(i_pt_bin), -1, result["x"], result["mean_guess"], result["mean_guess_err"], result["mean_guess_err_bootstrap"], result["mean_fit"], result["mean_fit_err"], result["std_guess"], result["std_fit"], reco_pt_i]
        df_inclusive.loc[i_pt_bin] = result_i

    else:
        result_i = [int(i_pt_bin), -1, 0, 0, 0, 0, 0, 0, 0, reco_pt_i]
        #df_inclusive.loc[i_pt_bin] = result_i

for i_eta_bin in range(num_eta_bins):
    print(f"\nEta bin {i_eta_bin}")
    eta_low = eta_edges[i_eta_bin]
    eta_high = eta_edges[i_eta_bin+1]
    df_eta_bin = df[df["eta_bin"] == i_eta_bin]

    x = df_eta_bin["mean_pt"].values
    x_reco = df_eta_bin["reco_pt"].values
    y_fit = df_eta_bin["mean_scale_fit"].values
    y_hist = df_eta_bin["mean_scale"].values
    #dy_fit = df_eta_bin["mean_scale_fit_err"].values
    dy_fit = df_eta_bin["mean_scale_err_bootstrap"].values
    dy_hist = df_eta_bin["mean_scale_err"].values
    std = df_eta_bin["std_scale"].values
    rms = df_eta_bin["rms_scale"].values

    ## Plotting
    fig, ax = plt.subplots()
    ax.errorbar(x, y_fit, yerr=dy_fit, fmt="o", color="red", label="Fit")
    ax.errorbar(x, y_hist, yerr=dy_hist, fmt="^", color="orange", label="Fit")
    ax.set_xlabel(r"Mean L1 $p_T$ [GeV]")
    ax.set_ylabel(r"Mean L1/Reco $p_T$")
    ax.axhline(1.0, color="black", linestyle="--")
    ax.text(0.65, 0.75, f"eta bin [{eta_edges[i_eta_bin]}, {eta_edges[i_eta_bin+1]}]", transform=ax.transAxes, fontsize=16*1.2)
    ax.set_xscale("log")
    maj = [100, 500]
    minor_labels = {20, 30, 50, 200, 300}
    ax.xaxis.set_major_locator(FixedLocator(maj))
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10) * 0.1))
    ax.xaxis.set_minor_formatter(FuncFormatter(minor_fmt))
    plt.tight_layout()
    ax.legend(fontsize=16*1.2)
    hep.cms.label(ax=ax, data=False, rlabel="Level-1 Scouting")
    plt.savefig(f"{output_dir}/jet_pt_scale_eta_{i_eta_bin}.png")

    # Wrt reco pT
    fig, ax = plt.subplots()
    ax.errorbar(x_reco, y_fit, yerr=dy_fit, fmt="o", color="red", label="Fit")
    ax.errorbar(x_reco, y_hist, yerr=dy_hist, fmt="^", color="orange", label="Hist")
    ax.set_xlabel(r"Mean L1 $p_T$ [GeV]")
    ax.set_ylabel(r"Mean L1/Reco $p_T$")
    ax.axhline(1.0, color="black", linestyle="--")
    ax.text(0.65, 0.75, f"eta bin [{eta_edges[i_eta_bin]}, {eta_edges[i_eta_bin+1]}]", transform=ax.transAxes, fontsize=16*1.2)
    ax.set_xscale("log")
    maj = [100, 500]
    minor_labels = {20, 30, 50, 200, 300}
    ax.xaxis.set_major_locator(FixedLocator(maj))
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10) * 0.1))
    ax.xaxis.set_minor_formatter(FuncFormatter(minor_fmt))
    plt.tight_layout()
    ax.legend(fontsize=16*1.2)
    hep.cms.label(ax=ax, data=False, rlabel="Level-1 Scouting")
    plt.savefig(f"{output_dir}/jet_reco_pt_scale_eta_{i_eta_bin}.png")

    ## Fitting
    #kernel_fit = (C(1.0) * Matern(length_scale=1.0, length_scale_bounds=(1e-3, 1e3), nu=1.5)) + WhiteKernel(noise_level=0.01, noise_level_bounds=(1e-5, 1e3)) 
    kernel_fit = (C(1.0) * Matern(length_scale=1.0, length_scale_bounds=(1e-3, 1e3), nu=1.5)) + WhiteKernel(noise_level=0.01, noise_level_bounds=(1e-5, 1e3)) 
    #kernel_fit = (C(1.0) * Matern(length_scale=1.0, length_scale_bounds=(1e-3, 1e3), nu=1.5))
    gpr_fit = GaussianProcessRegressor(kernel=kernel_fit, alpha=dy_fit**2)
    gpr_fit.fit(np.log(x.reshape(-1, 1)), y_fit.reshape(-1, 1))

    #kernel_hist = (C(1.0) * Matern(length_scale=1.0, length_scale_bounds=(1e-3, 1e3), nu=1.5)) + WhiteKernel(noise_level=0.01, noise_level_bounds=(1e-5, 1e3)) 
    kernel_hist = (C(1.0) * Matern(length_scale=1.0, length_scale_bounds=(1e-3, 1e3), nu=1.5)) + WhiteKernel(noise_level=0.01, noise_level_bounds=(1e-5, 1e3)) 
    #kernel_hist = (C(1.0) * Matern(length_scale=1.0, length_scale_bounds=(1e-3, 1e3), nu=1.5))
    gpr_hist = GaussianProcessRegressor(kernel=kernel_hist, alpha=dy_hist**2)
    gpr_hist.fit(np.log(x.reshape(-1, 1)), y_hist.reshape(-1, 1))

    x_plot = np.linspace(20, 600, 1000)
    x_plot_log = np.log(x_plot.reshape(-1, 1))
    y_fit_plot, dy_fit_plot = gpr_fit.predict(x_plot_log, return_std=True)
    y_hist_plot, dy_hist_plot = gpr_hist.predict(x_plot_log, return_std=True)

    # Plot the fits
    fig, ax = plt.subplots()
    ax.errorbar(x, y_fit, yerr=dy_fit, fmt="o", color="red", label="Fit")
    ax.plot(x_plot, y_fit_plot, '-', label="Fit GPR", color="red")
    ax.fill_between(x_plot, y_fit_plot - dy_fit_plot, y_fit_plot + dy_fit_plot, color='red', alpha=0.2)
    ax.errorbar(x, y_hist, yerr=dy_hist, fmt="^", color="orange", label="Hist")
    ax.plot(x_plot, y_hist_plot, '-', label="Hist GPR", color="orange")
    ax.fill_between(x_plot, y_hist_plot - dy_hist_plot, y_hist_plot + dy_hist_plot, color='orange', alpha=0.2)
    ax.set_xlabel(r"Mean L1 $p_T$ [GeV]")
    ax.set_ylabel(r"Mean L1/Reco $p_T$")
    ax.axhline(1.0, color="black", linestyle="--")
    ax.text(0.65, 0.75, f"eta bin [{eta_edges[i_eta_bin]}, {eta_edges[i_eta_bin+1]}]", transform=ax.transAxes, fontsize=16*1.2)
    ax.set_xscale("log")
    maj = [100, 500]
    minor_labels = {20, 30, 50, 200, 300}
    ax.xaxis.set_major_locator(FixedLocator(maj))
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10) * 0.1))
    ax.xaxis.set_minor_formatter(FuncFormatter(minor_fmt))
    plt.tight_layout()
    ax.legend(fontsize=16*1.2)
    hep.cms.label(ax=ax, data=False, rlabel="Level-1 Scouting")
    plt.savefig(f"{output_dir}/jet_pt_scale_fit_eta_{i_eta_bin}.png")
