# Extract the scale response from the fits, plot and fit a GPR to the scale for two datasets
# Also plots the ratio of the scale responses and saves the histograms
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
parser.add_argument("--input_a", "-ia", type=str, help="Input directory of first series of fits")
parser.add_argument("--label_a", "-la", type=str, help="Label of first series of fits")
parser.add_argument("--input_b", "-ib", type=str, help="Input directory of second series of fits")
parser.add_argument("--label_b", "-lb", type=str, help="Label of second series of fits")
parser.add_argument("--output", "-o", type=str, help="Output plot directory (also stores fits)")
args = parser.parse_args()

input_dir_a = args.input_a
label_a = args.label_a
input_dir_b = args.input_b
label_b = args.label_b
output_dir = args.output
cwd = os.getcwd()

os.makedirs(output_dir, exist_ok=True)

###### Histogram initializations ########
# Store hists that give the data/MC correction
eta_axis = hist.axis.Variable([-3.0, -2.5, -2.0, -1.3, -0.5, 0.0, 0.5, 1.3, 2.0, 2.5, 3.0], name="eta", label="Jet eta")
pt_axis = hist.axis.Regular(215, 20, 450, name="pt", label="Jet pT [GeV]")

# The individual pT corrections
h_input_a_correction_fit_nominal = hist.Hist(eta_axis, pt_axis, storage=hist.storage.Weight())
h_input_a_correction_fit_up = hist.Hist(eta_axis, pt_axis, storage=hist.storage.Weight())
h_input_a_correction_fit_down = hist.Hist(eta_axis, pt_axis, storage=hist.storage.Weight())

h_input_a_correction_hist_nominal = hist.Hist(eta_axis, pt_axis, storage=hist.storage.Weight())
h_input_a_correction_hist_up = hist.Hist(eta_axis, pt_axis, storage=hist.storage.Weight())
h_input_a_correction_hist_down = hist.Hist(eta_axis, pt_axis, storage=hist.storage.Weight())

h_input_b_correction_fit_nominal = hist.Hist(eta_axis, pt_axis, storage=hist.storage.Weight())
h_input_b_correction_fit_up = hist.Hist(eta_axis, pt_axis, storage=hist.storage.Weight())
h_input_b_correction_fit_down = hist.Hist(eta_axis, pt_axis, storage=hist.storage.Weight())

h_input_b_correction_hist_nominal = hist.Hist(eta_axis, pt_axis, storage=hist.storage.Weight())
h_input_b_correction_hist_up = hist.Hist(eta_axis, pt_axis, storage=hist.storage.Weight())
h_input_b_correction_hist_down = hist.Hist(eta_axis, pt_axis, storage=hist.storage.Weight())

# Scale corrections
h_scale_correction_fit_nominal = hist.Hist(eta_axis, pt_axis, storage=hist.storage.Weight())
h_scale_correction_fit_up = hist.Hist(eta_axis, pt_axis, storage=hist.storage.Weight())
h_scale_correction_fit_down = hist.Hist(eta_axis, pt_axis, storage=hist.storage.Weight())
scale_pt_centers = h_scale_correction_fit_nominal.axes[1].centers

h_scale_correction_hist_nominal = hist.Hist(eta_axis, pt_axis, storage=hist.storage.Weight())
h_scale_correction_hist_up = hist.Hist(eta_axis, pt_axis, storage=hist.storage.Weight())
h_scale_correction_hist_down = hist.Hist(eta_axis, pt_axis, storage=hist.storage.Weight())
scale_pt_centers = h_scale_correction_hist_nominal.axes[1].centers

###### Data filling ########
# Make sure these are consistent with production
pt_edges = np.concatenate([np.arange(20, 60, 4), np.arange(60, 100, 5), np.arange(100, 200, 10), np.arange(200, 300, 20), np.arange(300, 500, 50)])
eta_edges = np.array([-3.0, -2.5, -2.0, -1.3, -0.5, 0.0, 0.5, 1.3, 2.0, 2.5, 3.0])
num_pt_bins = len(pt_edges) - 1
num_eta_bins = len(eta_edges) - 1

# Make a pandas dataframe to collect the fit results to plot
col_names = ["pt_bin", "eta_bin", "mean_pt", "mean_scale", "mean_scale_err", "mean_scale_err_bootstrap", "mean_scale_fit", "mean_scale_fit_err", "rms_scale", "std_scale", "reco_pt"]
df_a = pd.DataFrame(columns=col_names)
df_b= pd.DataFrame(columns=col_names)

# Loop over all the eta bins first, then pT bins and fill dataframe
i_row = 0
for i_eta_bin in range(num_eta_bins):
    for i_pt_bin in range(num_pt_bins):
    #for i_pt_bin in range(2, num_pt_bins):
        fname_a = f"{input_dir_a}/fit_pt_{i_pt_bin}_eta_{i_eta_bin}.pkl"
        fname_b = f"{input_dir_b}/fit_pt_{i_pt_bin}_eta_{i_eta_bin}.pkl"

        reco_pt_i = 0.5 * (pt_edges[i_pt_bin] + pt_edges[i_pt_bin+1])

        result_a = None
        result_b = None

        # Load first sample
        if os.path.exists(fname_a):
            with open(fname_a, "rb") as f:
                result_a = pickle.load(f)
                result_a = result_a["scale"]
 
        if result_a is not None:
            result_a_i = [int(i_pt_bin), int(i_eta_bin), result_a["x"], result_a["mean_guess"], result_a["mean_guess_err"], result_a["mean_guess_err_bootstrap"], result_a["mean_fit"], result_a["mean_fit_err"], result_a["std_guess"], result_a["std_fit"], reco_pt_i]
            df_a.loc[i_row] = result_a_i

        else:
            result_a_i = [int(i_pt_bin), int(i_eta_bin), 0, 0, 0, 0, 0, 0, 0, reco_pt_i]
            #df.loc[i_row] = result_a_i

        # Load second sample
        if os.path.exists(fname_b):
            with open(fname_b, "rb") as f:
                result_b = pickle.load(f)
                result_b = result_b["scale"]
 
        if result_b is not None:
            result_b_i = [int(i_pt_bin), int(i_eta_bin), result_b["x"], result_b["mean_guess"], result_b["mean_guess_err"], result_b["mean_guess_err_bootstrap"], result_b["mean_fit"], result_b["mean_fit_err"], result_b["std_guess"], result_b["std_fit"], reco_pt_i]
            df_b.loc[i_row] = result_b_i

        else:
            result_b_i = [int(i_pt_bin), int(i_eta_bin), 0, 0, 0, 0, 0, 0, 0, reco_pt_i]
            #df.loc[i_row] = result_b_i
            
        i_row += 1

######### Plot, fit and extract ###############
for i_eta_bin in range(num_eta_bins):
    print(f"\nEta bin {i_eta_bin}")
    eta_low = eta_edges[i_eta_bin]
    eta_high = eta_edges[i_eta_bin+1]

    df_a_eta_bin = df_a[df_a["eta_bin"] == i_eta_bin]
    x_a = df_a_eta_bin["mean_pt"].values
    x_reco_a = df_a_eta_bin["reco_pt"].values
    y_fit_a = df_a_eta_bin["mean_scale_fit"].values
    y_hist_a = df_a_eta_bin["mean_scale"].values
    #dy_fit_a = df_a_eta_bin["mean_scale_fit_err"].values
    dy_fit_a = df_a_eta_bin["mean_scale_err_bootstrap"].values
    dy_hist_a = df_a_eta_bin["mean_scale_err"].values
    std_a = df_a_eta_bin["std_scale"].values
    rms_a = df_a_eta_bin["rms_scale"].values

    df_b_eta_bin = df_b[df_b["eta_bin"] == i_eta_bin]
    x_b = df_b_eta_bin["mean_pt"].values
    x_reco_b = df_b_eta_bin["reco_pt"].values
    y_fit_b = df_b_eta_bin["mean_scale_fit"].values
    y_hist_b = df_b_eta_bin["mean_scale"].values
    #dy_fit_b = df_b_eta_bin["mean_scale_fit_err"].values
    dy_fit_b = df_b_eta_bin["mean_scale_err_bootstrap"].values
    dy_hist_b = df_b_eta_bin["mean_scale_err"].values
    std_b = df_b_eta_bin["std_scale"].values
    rms_b = df_b_eta_bin["rms_scale"].values

    ## Plotting (fitted means)
    fig, ax = plt.subplots()
    ax.errorbar(x_a, y_fit_a, yerr=dy_fit_a, fmt="o", color="red", label=f"{label_a} Fit")
    ax.errorbar(x_a, y_fit_b, yerr=dy_fit_b, fmt="o", color="blue", label=f"{label_b} Fit")
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

    ## Plotting (histogram means)
    fig, ax = plt.subplots()
    ax.errorbar(x_a, y_hist_a, yerr=dy_hist_a, fmt="o", color="red", label=f"{label_a} Hist")
    ax.errorbar(x_a, y_hist_b, yerr=dy_hist_b, fmt="o", color="blue", label=f"{label_b} Hist")
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
    plt.savefig(f"{output_dir}/jet_pt_scale_hist_eta_{i_eta_bin}.png")
    
    ## Fitting
    # a
    kernel_fit_a = (C(1.0) * Matern(length_scale=1.0, length_scale_bounds=(1e-3, 1e3), nu=1.5)) + WhiteKernel(noise_level=0.01, noise_level_bounds=(1e-5, 1e3)) 
    gpr_fit_a = GaussianProcessRegressor(kernel=kernel_fit_a, alpha=dy_fit_a**2)
    gpr_fit_a.fit(np.log(x_a.reshape(-1, 1)), y_fit_a.reshape(-1, 1))

    kernel_hist_a = (C(1.0) * Matern(length_scale=1.0, length_scale_bounds=(1e-3, 1e3), nu=1.5)) + WhiteKernel(noise_level=0.01, noise_level_bounds=(1e-5, 1e3)) 
    gpr_hist_a = GaussianProcessRegressor(kernel=kernel_hist_a, alpha=dy_hist_a**2)
    gpr_hist_a.fit(np.log(x_a.reshape(-1, 1)), y_hist_a.reshape(-1, 1))

    x_a_plot = np.linspace(20, 600, 1000)
    x_a_plot_log = np.log(x_a_plot.reshape(-1, 1))
    y_fit_a_plot, dy_fit_a_plot = gpr_fit_a.predict(x_a_plot_log, return_std=True)
    y_hist_a_plot, dy_hist_a_plot = gpr_hist_a.predict(x_a_plot_log, return_std=True)

    # b
    kernel_fit_b = (C(1.0) * Matern(length_scale=1.0, length_scale_bounds=(1e-3, 1e3), nu=1.5)) + WhiteKernel(noise_level=0.01, noise_level_bounds=(1e-5, 1e3)) 
    gpr_fit_b = GaussianProcessRegressor(kernel=kernel_fit_b, alpha=dy_fit_b**2)
    gpr_fit_b.fit(np.log(x_b.reshape(-1, 1)), y_fit_b.reshape(-1, 1))

    kernel_hist_b = (C(1.0) * Matern(length_scale=1.0, length_scale_bounds=(1e-3, 1e3), nu=1.5)) + WhiteKernel(noise_level=0.01, noise_level_bounds=(1e-5, 1e3)) 
    gpr_hist_b = GaussianProcessRegressor(kernel=kernel_hist_b, alpha=dy_hist_b**2)
    gpr_hist_b.fit(np.log(x_b.reshape(-1, 1)), y_hist_b.reshape(-1, 1))

    x_b_plot = np.linspace(20, 600, 1000)
    x_b_plot_log = np.log(x_b_plot.reshape(-1, 1))
    y_fit_b_plot, dy_fit_b_plot = gpr_fit_b.predict(x_b_plot_log, return_std=True)
    y_hist_b_plot, dy_hist_b_plot = gpr_hist_b.predict(x_b_plot_log, return_std=True)

    ## Ratios 
    ratio_fit = y_fit_a_plot/y_fit_b_plot
    ratio_fit_uncert = ratio_fit * np.sqrt((dy_fit_a_plot/y_fit_a_plot)**2 + (dy_fit_b_plot/y_fit_b_plot)**2)

    ratio_hist = y_hist_a_plot/y_hist_b_plot
    ratio_hist_uncert = ratio_hist * np.sqrt((dy_hist_a_plot/y_hist_a_plot)**2 + (dy_hist_b_plot/y_hist_b_plot)**2)

    ## Plotting GPR (fitted means)
    fig, axs = plt.subplots(2, 1, gridspec_kw=dict(height_ratios=[3, 1], hspace=0.1), sharex=True)
    ax.errorbar(x_a, y_fit_a, yerr=dy_fit_a, fmt="o", color="red", label=f"{label_a} Fit")
    ax.plot(x_a_plot, y_fit_a_plot, '-', label=f"{label_a} Fit GPR", color="red")
    ax.fill_between(x_a_plot, y_fit_a_plot - dy_fit_a_plot, y_fit_a_plot + dy_fit_a_plot, color='red', alpha=0.2)
    ax.errorbar(x_b, y_fit_b, yerr=dy_fit_b, fmt="o", color="blue", label=f"{label_b} Fit")
    ax.plot(x_b_plot, y_fit_b_plot, '-', label=f"{label_b} Fit GPR", color="blue")
    ax.fill_between(x_b_plot, y_fit_b_plot - dy_fit_b_plot, y_fit_b_plot + dy_fit_b_plot, color='blue', alpha=0.2)
    axs[0].set_xlabel("")
    axs[0].set_ylabel(r"Mean L1/Reco $p_T$")
    axs[0].axhline(1.0, color="black", linestyle="--")
    axs[0].legend(fontsize=16*1.2)
    #axs[0].set_ylim(0.7, 1.3)
    axs[1].text(0.70, 0.70, f"eta bin [{eta_edges[i_eta_bin]}, {eta_edges[i_eta_bin+1]}]", transform=axs[1].transAxes, fontsize=16*1.2)
    hep.cms.label(ax=axs[0], data=False, rlabel="Level-1 Scouting")
    # Ratio
    axs[1].plot(x_a_plot, ratio_fit, "-", label="Data/MC", color="red")
    axs[1].fill_between(x_a_plot, ratio_fit - ratio_fit_uncert, ratio_fit + ratio_fit_uncert, color='red', alpha=0.2)
    axs[1].axhline(1.0, color="black", linestyle="--")
    axs[1].axvline(30, color="black", linestyle="--")
    axs[0].axvline(30, color="black", linestyle="--")
    axs[1].set_xlabel(r"Mean L1 $p_T$ [GeV]")
    axs[1].set_ylabel(f"{label_a}/{label_b}")
    axs[1].set_ylim(0.8, 1.2)
    axs[0].set_xscale("log")
    maj = [100, 500]
    minor_labels = {20, 30, 50, 200, 300}
    axs[1].xaxis.set_major_locator(FixedLocator(maj))
    axs[1].xaxis.set_major_formatter(ScalarFormatter())
    axs[1].xaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10) * 0.1))
    axs[1].xaxis.set_minor_formatter(FuncFormatter(minor_fmt))
    plt.tight_layout()
    plt.savefig(f"{output_dir}/jet_pt_scale_fit_gpr_ratio_eta_{i_eta_bin}.png")

    ## Plotting GPR (histogram means)
    fig, axs = plt.subplots(2, 1, gridspec_kw=dict(height_ratios=[3, 1], hspace=0.1), sharex=True)
    ax.errorbar(x_a, y_hist_a, yerr=dy_hist_a, fmt="o", color="red", label=f"{label_a} Hist")
    ax.plot(x_a_plot, y_hist_a_plot, '-', label=f"{label_a} Hist GPR", color="red")
    ax.fill_between(x_a_plot, y_hist_a_plot - dy_hist_a_plot, y_hist_a_plot + dy_hist_a_plot, color='red', alpha=0.2)
    ax.errorbar(x_b, y_hist_b, yerr=dy_hist_b, fmt="o", color="blue", label=f"{label_b} Hist")
    ax.plot(x_b_plot, y_hist_b_plot, '-', label=f"{label_b} Hist GPR", color="blue")
    ax.fill_between(x_b_plot, y_hist_b_plot - dy_hist_b_plot, y_hist_b_plot + dy_hist_b_plot, color='blue', alpha=0.2)
    axs[0].set_xlabel("")
    axs[0].set_ylabel(r"Mean L1/Reco $p_T$")
    axs[0].axhline(1.0, color="black", linestyle="--")
    axs[0].legend(fontsize=16*1.2)
    #axs[0].set_ylim(0.7, 1.3)
    axs[1].text(0.70, 0.70, f"eta bin [{eta_edges[i_eta_bin]}, {eta_edges[i_eta_bin+1]}]", transform=axs[1].transAxes, fontsize=16*1.2)
    hep.cms.label(ax=axs[0], data=False, rlabel="Level-1 Scouting")
    # Ratio
    axs[1].plot(x_a_plot, ratio_hist, "-", label="Data/MC", color="red")
    axs[1].fill_between(x_a_plot, ratio_hist - ratio_hist_uncert, ratio_hist + ratio_hist_uncert, color='red', alpha=0.2)
    axs[1].axhline(1.0, color="black", linestyle="--")
    axs[1].axvline(30, color="black", linestyle="--")
    axs[0].axvline(30, color="black", linestyle="--")
    axs[1].set_xlabel(r"Mean L1 $p_T$ [GeV]")
    axs[1].set_ylabel(f"{label_a}/{label_b}")
    axs[1].set_ylim(0.8, 1.2)
    axs[0].set_xscale("log")
    maj = [100, 500]
    minor_labels = {20, 30, 50, 200, 300}
    axs[1].xaxis.set_major_locator(FixedLocator(maj))
    axs[1].xaxis.set_major_formatter(ScalarFormatter())
    axs[1].xaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10) * 0.1))
    axs[1].xaxis.set_minor_formatter(FuncFormatter(minor_fmt))
    plt.tight_layout()
    plt.savefig(f"{output_dir}/jet_pt_scale_hist_gpr_ratio_eta_{i_eta_bin}.png")

    ## Save the histograms
    scale_fit_a, scale_fit_err_a = gpr_fit_a.predict(np.log(scale_pt_centers.reshape(-1, 1)), return_std=True)
    scale_hist_a, scale_hist_err_a = gpr_hist_a.predict(np.log(scale_pt_centers.reshape(-1, 1)), return_std=True)

    scale_fit_b, scale_fit_err_b = gpr_fit_b.predict(np.log(scale_pt_centers.reshape(-1, 1)), return_std=True)
    scale_hist_b, scale_hist_err_b = gpr_hist_b.predict(np.log(scale_pt_centers.reshape(-1, 1)), return_std=True)

    scale_correction_fit = scale_fit_a / scale_fit_b
    scale_correction_fit_err = scale_correction_fit * np.sqrt((scale_fit_err_a/scale_fit_a)**2 + (scale_fit_err_b/scale_fit_b)**2)
    scale_correction_fit_up = scale_correction_fit + scale_correction_fit_err
    scale_correction_fit_down = scale_correction_fit - scale_correction_fit_err
    scale_correction_fit = scale_correction_fit.ravel()
    scale_correction_fit_up = scale_correction_fit_up.ravel()
    scale_correction_fit_down = scale_correction_fit_down.ravel()

    scale_correction_hist = scale_hist_a / scale_hist_b
    scale_correction_hist_err = scale_correction_hist * np.sqrt((scale_hist_err_a/scale_hist_a)**2 + (scale_hist_err_b/scale_hist_b)**2)
    scale_correction_hist_up = scale_correction_hist + scale_correction_hist_err
    scale_correction_hist_down = scale_correction_hist - scale_correction_hist_err
    scale_correction_hist = scale_correction_hist.ravel()
    scale_correction_hist_up = scale_correction_hist_up.ravel()
    scale_correction_hist_down = scale_correction_hist_down.ravel()

    ## a
    # Fit
    scale_fit_a = scale_fit_a.ravel()
    scale_fit_err_a = scale_fit_err_a.ravel()
    scale_fit_up_a = scale_fit_a + scale_fit_err_a
    scale_fit_down_a = scale_fit_a - scale_fit_err_a
    h_input_a_correction_fit_nominal.view().value[i_eta_bin, :] = scale_fit_a
    h_input_a_correction_fit_up.view().value[i_eta_bin, :] = scale_fit_up_a
    h_input_a_correction_fit_down.view().value[i_eta_bin, :] = scale_fit_down_a

    # Hist
    scale_hist_a = scale_hist_a.ravel()
    scale_hist_err_a = scale_hist_err_a.ravel()
    scale_hist_up_a = scale_hist_a + scale_hist_err_a
    scale_hist_down_a = scale_hist_a - scale_hist_err_a
    h_input_a_correction_hist_nominal.view().value[i_eta_bin, :] = scale_hist_a
    h_input_a_correction_hist_up.view().value[i_eta_bin, :] = scale_hist_up_a
    h_input_a_correction_hist_down.view().value[i_eta_bin, :] = scale_hist_down_a

    ## b
    # Fit
    scale_fit_b = scale_fit_b.ravel()
    scale_fit_err_b = scale_fit_err_b.ravel()
    scale_fit_up_b = scale_fit_b + scale_fit_err_b
    scale_fit_down_b = scale_fit_b - scale_fit_err_b
    h_input_b_correction_fit_nominal.view().value[i_eta_bin, :] = scale_fit_b
    h_input_b_correction_fit_up.view().value[i_eta_bin, :] = scale_fit_up_b
    h_input_b_correction_fit_down.view().value[i_eta_bin, :] = scale_fit_down_b

    # Hist
    scale_hist_b = scale_hist_b.ravel()
    scale_hist_err_b = scale_hist_err_b.ravel()
    scale_hist_up_b = scale_hist_b + scale_hist_err_b
    scale_hist_down_b = scale_hist_b - scale_hist_err_b
    h_input_b_correction_hist_nominal.view().value[i_eta_bin, :] = scale_hist_b
    h_input_b_correction_hist_up.view().value[i_eta_bin, :] = scale_hist_up_b
    h_input_b_correction_hist_down.view().value[i_eta_bin, :] = scale_hist_down_b

    ## Ratio
    # Fit
    h_scale_correction_fit_nominal.view().value[i_eta_bin, :] = scale_correction_fit
    h_scale_correction_fit_up.view().value[i_eta_bin, :] = scale_correction_fit_up
    h_scale_correction_fit_down.view().value[i_eta_bin, :] = scale_correction_fit_down

    # Hist
    h_scale_correction_hist_nominal.view().value[i_eta_bin, :] = scale_correction_hist
    h_scale_correction_hist_up.view().value[i_eta_bin, :] = scale_correction_hist_up
    h_scale_correction_hist_down.view().value[i_eta_bin, :] = scale_correction_hist_down

# Save the histograms 
with uproot.recreate(f"{output_dir}/scale_{label_a}_fit.root") as f:
    f["h_scale_nominal"] = h_input_a_correction_fit_nominal
    f["h_scale_up"] = h_input_a_correction_fit_up
    f["h_scale_down"] = h_input_a_correction_fit_down

with uproot.recreate(f"{output_dir}/scale_{label_a}_hist.root") as f:
    f["h_scale_nominal"] = h_input_a_correction_hist_nominal
    f["h_scale_up"] = h_input_a_correction_hist_up
    f["h_scale_down"] = h_input_a_correction_hist_down

with uproot.recreate(f"{output_dir}/scale_{label_b}_fit.root") as f:
    f["h_scale_nominal"] = h_input_b_correction_fit_nominal
    f["h_scale_up"] = h_input_b_correction_fit_up
    f["h_scale_down"] = h_input_b_correction_fit_down

with uproot.recreate(f"{output_dir}/scale_{label_b}_hist.root") as f:
    f["h_scale_nominal"] = h_input_b_correction_hist_nominal
    f["h_scale_up"] = h_input_b_correction_hist_up
    f["h_scale_down"] = h_input_b_correction_hist_down

with uproot.recreate(f"{output_dir}/scale_{label_a}_{label_b}_ratio_fit.root") as f:
    f["h_scale_nominal"] = h_scale_correction_fit_nominal
    f["h_scale_up"] = h_scale_correction_fit_up
    f["h_scale_down"] = h_scale_correction_fit_down

with uproot.recreate(f"{output_dir}/scale_{label_a}_{label_b}_ratio_hist.root") as f:
    f["h_scale_nominal"] = h_scale_correction_hist_nominal
    f["h_scale_up"] = h_scale_correction_hist_up
    f["h_scale_down"] = h_scale_correction_hist_down