# This script extracts the mean pT, and quantities required for the scale and resolution
# Example command: python3 Preprocess/fit_histo.py --input outputs/Histograms/histos_eta_inclusive_pt_3.root --output outputs/Fits/test/results.pkl --ptlow 50 --pthigh 60 --etalow neg3p0 --etahigh pos3p0
import os
import argparse
import ROOT

import uproot
import hist
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit

import matplotlib.pyplot as plt
import mplhep as hep
hep.style.use("CMS")

# Encoding function to convert -3.0 to "neg3p0" and 3.0 to "pos3p0"
def encode_eta(eta):
    if eta < 0:
        return "neg" + str(abs(eta)).replace(".", "p")
    else:
        return "pos" + str(eta).replace(".", "p")

# A converter to decode "pos3p0" to 3.0 and "neg3p0" to -3.0
def decode_eta(eta_str):
    if eta_str.startswith("pos"):
        return float(eta_str[3:].replace("p", "."))
    elif eta_str.startswith("neg"):
        return -float(eta_str[3:].replace("p", "."))
    else:
        return eta_str

#Argument parser
parser = argparse.ArgumentParser(description='Create histograms')
parser.add_argument('--input', '-i', type=str, help='Input file')
parser.add_argument('--output', '-o', type=str, help='Output file')
parser.add_argument('--ptlow', type=str, help='pT bin low edge')
parser.add_argument('--pthigh', type=str, help='pT bin high edge')
parser.add_argument('--etalow', type=str, help='eta bin low edge')
parser.add_argument('--etahigh', type=str, help='eta bin high edge')
parser.add_argument('--rebin', type=int, help="Rebin factor")
parser.add_argument('--plotdir', '-p', type=str, default="plots/Fit/data", help='Plot directory')
args = parser.parse_args()

infile = args.input
outfile = args.output
rebin = args.rebin
# Metadata
ptlow = float(args.ptlow)
pthigh = float(args.pthigh)
etalow = args.etalow
etahigh = args.etahigh
# Decode
etalow = decode_eta(etalow)
etahigh = decode_eta(etahigh)


# Gaussian fit function
def gauss(x, a, mean, sig):
    return a * np.exp(-(x - mean)**2 / (2 * sig**2))

f = uproot.open(infile)
h_pt_scale = f["h_pt_scale"].to_hist()
h_pt_reso = f["h_pt_reso"].to_hist()
h_diff = f["h_diff"].to_hist()
h_scale = f["h_scale"].to_hist()

rebin_factor = 1
if rebin is not None: 
    #h_diff = h_diff[::(rebin*1j)]
    h_scale = h_scale[::(rebin*1j)]
    rebin_factor = rebin
    
# Get the mean pT for the scale and resolution calculation
pt_scale = h_pt_scale.axes[0].centers
pt_scale_values = h_pt_scale.values()
pt_reso = h_pt_reso.axes[0].centers
pt_reso_values = h_pt_reso.values()
mean_pt_scale = np.average(pt_scale, weights=pt_scale_values)
mean_pt_reso = np.average(pt_reso, weights=pt_reso_values)

x_diff = h_diff.axes[0].centers
x_scale = h_scale.axes[0].centers
y_diff = h_diff.values()
y_scale = h_scale.values()
dy_diff = np.sqrt(y_diff)
dy_scale = np.sqrt(y_scale)

####### Scale fit #########
center_pos_scale = np.argmax(y_scale)
x_scale_pos_low = np.max([center_pos_scale - (6 // rebin_factor), 0])
x_scale_pos_high = np.min([center_pos_scale + (6 // rebin_factor), len(x_scale)])
fit_x_scale = x_scale[x_scale_pos_low:x_scale_pos_high]
fit_y_scale = y_scale[x_scale_pos_low:x_scale_pos_high]
fit_dy_scale = dy_scale[x_scale_pos_low:x_scale_pos_high]
# Guess parameters
a_scale_guess = y_scale[center_pos_scale]
mean_scale_guess = x_scale[center_pos_scale]
sig_scale_guess = np.average((fit_x_scale - mean_scale_guess)**2, weights=fit_y_scale)**0.5

#print(fit_x_scale)
#print(fit_y_scale)
#print(f"Guess parameters for Gaussian fit to scale: a = {a_scale_guess:.2f}, mean = {mean_scale_guess:.2f}, sig = {sig_scale_guess:.2f}")

# Fit to Gaussian
try:
    p0_scale = [a_scale_guess, mean_scale_guess, sig_scale_guess]
    popt_scale, pcov_scale = curve_fit(gauss, fit_x_scale, fit_y_scale, p0=p0_scale, sigma=fit_dy_scale, absolute_sigma=True)
    print(popt_scale)
    print(pcov_scale)

    # Plot the fit
    plot_x_scale_gauss = np.linspace(fit_x_scale[0], fit_x_scale[-1], 100)
    plot_y_scale_gauss = gauss(plot_x_scale_gauss, *popt_scale)
    fig, ax = plt.subplots()
    hep.histplot(h_scale, ax=ax, label="Data", color="black", histtype="step", flow=None)
    ax.plot(plot_x_scale_gauss, plot_y_scale_gauss, label="Gaussian Fit", color="red")
    ax.set_xlabel("L1/RecoJet pT")
    ax.set_ylabel("Events")
    hep.cms.label(ax=ax, data=False, rlabel="Level-1 Scouting 2024")
    ax.legend()
    ax.text(0.65, 0.75, f"pT bin [{ptlow}, {pthigh}]", transform=ax.transAxes, fontsize=16*1.2)
    ax.text(0.65, 0.70, f"eta bin [{etalow}, {etahigh}]", transform=ax.transAxes, fontsize=16*1.2)
    plt.savefig(f"{args.plotdir}/fit_scale_pt_{int(ptlow)}_{int(pthigh)}_eta_{encode_eta(etalow)}_{encode_eta(etahigh)}.png")

    ## Histogram stats
    # Get mean and rms from the whole region
    mean_scale = np.average(x_scale, weights=y_scale)
    rms_scale = np.average((x_scale - mean_scale)**2, weights=y_scale)**0.5
    # Get mean and rms from the fitting region
    #mean_scale = np.average(fit_x_scale, weights=fit_y_scale)
    #rms_scale = np.average((fit_x_scale - mean_scale)**2, weights=fit_y_scale)**0.5

    # Using 68% CL interval
    #cumsum_scale = np.cumsum(y_scale)
    #cumsum_scale = cumsum_scale / cumsum_scale[-1]
    #cumsum_func = sp.interpolate.interp1d(x_scale, cumsum_scale)
    #x_16 = sp.optimize.brentq(lambda x: cumsum_func(x) - 0.16, x_scale[0], x_scale[-1])
    #x_84 = sp.optimize.brentq(lambda x: cumsum_func(x) - 0.84, x_scale[0], x_scale[-1])
    #mean_scale = sp.optimize.brentq(lambda x: cumsum_func(x) - 0.50, x_scale[0], x_scale[-1])
    #rms_scale = (x_84 - x_16) / 2.0

    # Bootstrap to get the errors on the mean and RMS with fitting
    mean_scale_bootstrap = []
    for i_bootstrap in range(100):
        y_scale_bootstrap = np.random.poisson(y_scale)

        x_scale_pos_bootstrap = np.argmax(y_scale_bootstrap)
        x_scale_pos_low_bootstrap = np.max([x_scale_pos_bootstrap - (6 // rebin_factor), 0])
        x_scale_pos_high_bootstrap = np.min([x_scale_pos_bootstrap + (6 // rebin_factor), len(x_scale)])
        fit_x_scale_bootstrap = x_scale[x_scale_pos_low_bootstrap:x_scale_pos_high_bootstrap]
        fit_y_scale_bootstrap = y_scale_bootstrap[x_scale_pos_low_bootstrap:x_scale_pos_high_bootstrap]
        fit_dy_scale_bootstrap = dy_scale[x_scale_pos_low_bootstrap:x_scale_pos_high_bootstrap]

        #mean_scale_bootstrap.append(np.average(fit_x_scale_bootstrap, weights=fit_y_scale_bootstrap)) # No fit
        a_scale_guess_bootstrap = y_scale_bootstrap[x_scale_pos_bootstrap]
        mean_scale_guess_bootstrap = x_scale[x_scale_pos_bootstrap]
        sig_scale_guess_bootstrap = np.average((fit_x_scale_bootstrap - mean_scale_guess_bootstrap)**2, weights=fit_y_scale_bootstrap)**0.5
        try:
            p0_scale_bootstrap = [a_scale_guess_bootstrap, mean_scale_guess_bootstrap, sig_scale_guess_bootstrap]
            popt_scale_bootstrap, pcov_scale_bootstrap = curve_fit(gauss, fit_x_scale_bootstrap, fit_y_scale_bootstrap, p0=p0_scale_bootstrap, sigma=fit_dy_scale_bootstrap, absolute_sigma=True)
            mean_scale_bootstrap.append(popt_scale_bootstrap[1])
        except:
            print(f"Bootstrap fit {i_bootstrap} failed, skipping")

    mean_scale_bootstrap = np.array(mean_scale_bootstrap)
    print(f"Number of successful bootstraps: {len(mean_scale_bootstrap)}/100")

    mean_scale_fit = popt_scale[1]
    std_scale_fit = popt_scale[2]
    mean_scale_fit_err = np.sqrt(pcov_scale[1, 1])
    std_scale_fit_err = np.sqrt(pcov_scale[2, 2])
    # The error on the mean of the histogram from bootstrapping
    mean_scale_err_bootstrap = np.std(mean_scale_bootstrap)
    # The error on the mean of the histogram is just the RMS divided by sqrt(sum(y_scale))
    mean_scale_err = rms_scale/np.sqrt(np.sum(y_scale))
except:
    # Plot the histogram without fit
    fig, ax = plt.subplots()
    hep.histplot(h_scale, ax=ax, label="Data", color="black", histtype="step", flow=None)
    ax.set_xlabel("L1/RecoJet pT")
    ax.set_ylabel("Events")
    hep.cms.label(ax=ax, data=False, rlabel="Level-1 Scouting 2024")
    ax.legend()
    ax.text(0.65, 0.75, f"pT bin [{ptlow}, {pthigh}]", transform=ax.transAxes, fontsize=16*1.2)
    ax.text(0.65, 0.70, f"eta bin [{etalow}, {etahigh}]", transform=ax.transAxes, fontsize=16*1.2)
    plt.savefig(f"{args.plotdir}/fit_scale_pt_{int(ptlow)}_{int(pthigh)}_eta_{encode_eta(etalow)}_{encode_eta(etahigh)}.png")

    print(f"Scale fit failed for pt bin [{ptlow}, {pthigh}] and eta bin [{etalow}, {etahigh}]")

    ## Histogram stats
    # Get mean and rms from the whole region
    mean_scale = np.average(x_scale, weights=y_scale)
    rms_scale = np.average((x_scale - mean_scale)**2, weights=y_scale)**0.5
    # Get mean and rms from the fitting region
    #mean_scale = np.average(fit_x_scale, weights=fit_y_scale)
    #rms_scale = np.average((fit_x_scale - mean_scale)**2, weights=fit_y_scale)**0.5

    # Using 68% CL interval
    #cumsum_scale = np.cumsum(y_scale)
    #cumsum_scale = cumsum_scale / cumsum_scale[-1]
    #cumsum_func = sp.interpolate.interp1d(x_scale, cumsum_scale)
    #x_16 = sp.optimize.brentq(lambda x: cumsum_func(x) - 0.16, x_scale[0], x_scale[-1])
    #x_84 = sp.optimize.brentq(lambda x: cumsum_func(x) - 0.84, x_scale[0], x_scale[-1])
    #mean_scale = sp.optimize.brentq(lambda x: cumsum_func(x) - 0.50, x_scale[0], x_scale[-1])
    #rms_scale = (x_84 - x_16) / 2.0

    # Bootstrap to get the errors on the mean and RMS without fitting
    mean_scale_bootstrap = []
    for i_bootstrap in range(100):
        y_scale_bootstrap = np.random.poisson(y_scale)

        x_scale_pos_bootstrap = np.argmax(y_scale_bootstrap)
        x_scale_pos_low_bootstrap = np.max([x_scale_pos_bootstrap - (6 // rebin_factor), 0])
        x_scale_pos_high_bootstrap = np.min([x_scale_pos_bootstrap + (6 // rebin_factor), len(x_scale)])
        fit_x_scale_bootstrap = x_scale[x_scale_pos_low_bootstrap:x_scale_pos_high_bootstrap]
        fit_y_scale_bootstrap = y_scale_bootstrap[x_scale_pos_low_bootstrap:x_scale_pos_high_bootstrap]

        mean_scale_bootstrap.append(np.average(fit_x_scale_bootstrap, weights=fit_y_scale_bootstrap))

    mean_scale_bootstrap = np.array(mean_scale_bootstrap)
    print(f"Number of successful bootstraps: {len(mean_scale_bootstrap)}/100")

    mean_scale_fit = mean_scale
    std_scale_fit = rms_scale
    mean_scale_fit_err = rms_scale/np.sqrt(np.sum(y_scale))
    std_scale_fit_err = 0.0
    # The error on the mean of the histogram from bootstrapping
    mean_scale_err_bootstrap = np.std(mean_scale_bootstrap)
    # The error on the mean of the histogram is just the RMS divided by sqrt(sum(y_scale)) (TO DO: Compare with bootstrap)
    mean_scale_err = rms_scale/np.sqrt(np.sum(y_scale))
    

####### Resolution fit #########
center_pos_diff = np.argmax(y_diff)
x_diff_pos_low = np.max([center_pos_diff-5, 0])
x_diff_pos_high = np.min([center_pos_diff+6, len(x_diff)])
fit_x_diff = x_diff[x_diff_pos_low:x_diff_pos_high]
fit_y_diff = y_diff[x_diff_pos_low:x_diff_pos_high]
fit_dy_diff = dy_diff[x_diff_pos_low:x_diff_pos_high]
# Guess parameters
a_diff_guess = y_diff[center_pos_diff]
mean_diff_guess = x_diff[center_pos_diff]
sig_diff_guess = np.average((fit_x_diff - mean_diff_guess)**2, weights=fit_y_diff)**0.5

# Fit to Gaussian
try:
    p0_diff = [a_diff_guess, mean_diff_guess, sig_diff_guess]
    popt_diff, pcov_diff = curve_fit(gauss, fit_x_diff, fit_y_diff, p0=p0_diff, sigma=fit_dy_diff, absolute_sigma=True)
    print(popt_diff)
    print(pcov_diff)

    # Plot the fit
    plot_x_diff_gauss = np.linspace(fit_x_diff[0], fit_x_diff[-1], 100)
    plot_y_diff_gauss = gauss(plot_x_diff_gauss, *popt_diff)
    fig, ax = plt.subplots()
    hep.histplot(h_diff, ax=ax, label="Data", color="black", histtype="step", flow=None)
    ax.plot(plot_x_diff_gauss, plot_y_diff_gauss, label="Gaussian Fit", color="red")
    ax.set_xlabel("L1 - RecoJet pT (GeV)")
    ax.set_ylabel("Events")
    hep.cms.label(ax=ax, data=False, rlabel="Level-1 Scouting 2024")
    ax.legend()
    ax.text(0.65, 0.75, f"pT bin [{ptlow}, {pthigh}]", transform=ax.transAxes, fontsize=16*1.2)
    ax.text(0.65, 0.70, f"eta bin [{etalow}, {etahigh}]", transform=ax.transAxes, fontsize=16*1.2)
    plt.savefig(f"{args.plotdir}/fit_diff_pt_{int(ptlow)}_{int(pthigh)}_eta_{encode_eta(etalow)}_{encode_eta(etahigh)}.png")

    ## Histogram stats
    # Get mean and rms from the whole region
    #mean_diff = np.average(x_diff, weights=y_diff)
    #rms_diff = np.average((x_diff - mean_diff)**2, weights=y_diff)**0.5
    # Get mean and rms from the fitting region
    #mean_diff = np.average(fit_x_diff, weights=fit_y_diff)
    #rms_diff = np.average((fit_x_diff - mean_diff)**2, weights=fit_y_diff)**0.5

    # Using 68% CL interval
    cumsum_diff = np.cumsum(y_diff)
    cumsum_diff = cumsum_diff / cumsum_diff[-1]
    cumsum_func = sp.interpolate.interp1d(x_diff, cumsum_diff)
    x_16 = sp.optimize.brentq(lambda x: cumsum_func(x) - 0.16, x_diff[0], x_diff[-1])
    x_84 = sp.optimize.brentq(lambda x: cumsum_func(x) - 0.84, x_diff[0], x_diff[-1])
    mean_diff = sp.optimize.brentq(lambda x: cumsum_func(x) - 0.50, x_diff[0], x_diff[-1])
    #rms_diff = (x_84 - x_16) / 2.0
    # Take left sided interval to avoid falloff from the right
    rms_diff = mean_diff - x_16

    # Boostrap to get the error on the interval
    mean_diff_bootstrap = []
    rms_diff_bootstrap = []
    for i_bootstrap in range(100):
        y_diff_bootstrap = np.random.poisson(y_diff)
        cumsum_diff_bootstrap = np.cumsum(y_diff_bootstrap)
        cumsum_diff_bootstrap = cumsum_diff_bootstrap / cumsum_diff_bootstrap[-1]
        cumsum_func_bootstrap = sp.interpolate.interp1d(x_diff, cumsum_diff_bootstrap)
        x_16_bootstrap = sp.optimize.brentq(lambda x: cumsum_func_bootstrap(x) - 0.16, x_diff[0], x_diff[-1])
        x_84_bootstrap = sp.optimize.brentq(lambda x: cumsum_func_bootstrap(x) - 0.84, x_diff[0], x_diff[-1])
        x_50_diff_bootstrap = sp.optimize.brentq(lambda x: cumsum_func_bootstrap(x) - 0.50, x_diff[0], x_diff[-1])
        rms_diff_bootstrap.append(x_50_diff_bootstrap - x_16_bootstrap)
        mean_diff_bootstrap.append(x_50_diff_bootstrap)
    
    mean_diff_boostrap = np.array(mean_diff_bootstrap)
    rms_diff_bootstrap = np.array(rms_diff_bootstrap)


    mean_diff_fit = popt_diff[1]
    std_diff_fit = popt_diff[2]
    mean_diff_fit_err = np.sqrt(pcov_diff[1, 1])
    std_diff_fit_err = np.sqrt(pcov_diff[2, 2])
    # From bootstrapping
    mean_diff_err = np.std(mean_diff_bootstrap)
    rms_diff_err = np.std(rms_diff_bootstrap)

except:
    print(f"Resolution fit failed for pt bin [{ptlow}, {pthigh}] and eta bin [{etalow}, {etahigh}]")
    # Plot the histogram without fit
    fig, ax = plt.subplots()
    hep.histplot(h_diff, ax=ax, label="Data", color="black", histtype="step", flow=None)
    ax.set_xlabel("L1 - RecoJet pT (GeV)")
    ax.set_ylabel("Events")
    hep.cms.label(ax=ax, data=False, rlabel="Level-1 Scouting 2024")
    ax.legend()
    ax.text(0.65, 0.75, f"pT bin [{ptlow}, {pthigh}]", transform=ax.transAxes, fontsize=16*1.2)
    ax.text(0.65, 0.70, f"eta bin [{etalow}, {etahigh}]", transform=ax.transAxes, fontsize=16*1.2)
    plt.savefig(f"{args.plotdir}/fit_diff_pt_{int(ptlow)}_{int(pthigh)}_eta_{encode_eta(etalow)}_{encode_eta(etahigh)}.png")

    ## Histogram stats
    # Get mean and rms from the whole region
    #mean_diff = np.average(x_diff, weights=y_diff)
    #rms_diff = np.average((x_diff - mean_diff)**2, weights=y_diff)**0.5
    # Get mean and rms from the fitting region
    #mean_diff = np.average(fit_x_diff, weights=fit_y_diff)
    #rms_diff = np.average((fit_x_diff - mean_diff)**2, weights=fit_y_diff)**0.5

    # Using 68% CL interval
    cumsum_diff = np.cumsum(y_diff)
    cumsum_diff = cumsum_diff / cumsum_diff[-1]
    cumsum_func = sp.interpolate.interp1d(x_diff, cumsum_diff)
    x_16 = sp.optimize.brentq(lambda x: cumsum_func(x) - 0.16, x_diff[0], x_diff[-1])
    x_84 = sp.optimize.brentq(lambda x: cumsum_func(x) - 0.84, x_diff[0], x_diff[-1])
    mean_diff = sp.optimize.brentq(lambda x: cumsum_func(x) - 0.50, x_diff[0], x_diff[-1])
    #rms_diff = (x_84 - x_16) / 2.0
    # Take left sided interval to avoid falloff from the right
    rms_diff = mean_diff - x_16

    # Boostrap to get the error on the interval
    mean_diff_bootstrap = []
    rms_diff_bootstrap = []
    for i_bootstrap in range(100):
        y_diff_bootstrap = np.random.poisson(y_diff)
        cumsum_diff_bootstrap = np.cumsum(y_diff_bootstrap)
        cumsum_diff_bootstrap = cumsum_diff_bootstrap / cumsum_diff_bootstrap[-1]
        cumsum_func_bootstrap = sp.interpolate.interp1d(x_diff, cumsum_diff_bootstrap)
        x_16_bootstrap = sp.optimize.brentq(lambda x: cumsum_func_bootstrap(x) - 0.16, x_diff[0], x_diff[-1])
        x_84_bootstrap = sp.optimize.brentq(lambda x: cumsum_func_bootstrap(x) - 0.84, x_diff[0], x_diff[-1])
        x_50_diff_bootstrap = sp.optimize.brentq(lambda x: cumsum_func_bootstrap(x) - 0.50, x_diff[0], x_diff[-1])
        rms_diff_bootstrap.append(x_50_diff_bootstrap - x_16_bootstrap)
        mean_diff_bootstrap.append(x_50_diff_bootstrap)

    mean_diff_boostrap = np.array(mean_diff_bootstrap)
    rms_diff_bootstrap = np.array(rms_diff_bootstrap)

    mean_diff_fit = mean_diff
    std_diff_fit = rms_diff
    mean_diff_fit_err = np.std(mean_diff_bootstrap)
    std_diff_fit_err = np.std(rms_diff_bootstrap)
    # From bootstrapping
    mean_diff_err = np.std(mean_diff_bootstrap)
    rms_diff_err = np.std(rms_diff_bootstrap)

####### Save fit #########
print("\nFit results")
print("Scale")
print(f"  Mean x: {mean_pt_scale:.2f}")
print(f"  Guess Mean: {mean_scale:.2f}")
print(f"  Guess RMS: {rms_scale:.2f}")
print(f"  Fit Mean: {mean_scale_fit:.2f} +/- {mean_scale_fit_err:.2e}")
print(f"  Fit Std: {std_scale_fit:.2f} +/- {std_scale_fit_err:.2e}")
print("Resolution")
print(f"  Mean x: {mean_pt_reso:.2f}")
print(f"  Guess Mean: {mean_diff:.2f}")
print(f"  Guess RMS: {rms_diff:.2f}")
print(f"  Fit Mean: {mean_diff_fit:.2f} +/- {mean_diff_fit_err:.2e}")
print(f"  Fit Std: {std_diff_fit:.2f} +/- {std_diff_fit_err:.2e}")

import pickle
result_dict = {}
scale_dict = {}
resolution_dict = {}

scale_dict["x"] = mean_pt_scale
scale_dict["mean_guess"] = mean_scale
scale_dict["mean_guess_err"] = mean_scale_err
scale_dict["mean_guess_err_bootstrap"] = mean_scale_err_bootstrap
scale_dict["mean_fit"] = mean_scale_fit
scale_dict["mean_fit_err"] = mean_scale_fit_err
scale_dict["std_guess"] = rms_scale
scale_dict["std_fit"] = std_scale_fit
scale_dict["std_fit_err"] = std_scale_fit_err 

resolution_dict["x"] = mean_pt_reso
resolution_dict["mean_guess"] = mean_diff
resolution_dict["mean_guess_err"] = mean_diff_err
resolution_dict["mean_fit"] = mean_diff_fit
resolution_dict["mean_fit_err"] = mean_diff_fit_err
resolution_dict["std_guess"] = rms_diff
resolution_dict["std_guess_err"] = rms_diff_err
resolution_dict["std_fit"] = std_diff_fit
resolution_dict["std_fit_err"] = std_diff_fit_err 

result_dict["scale"] = scale_dict
result_dict["resolution"] = resolution_dict
# Add metadata
result_dict["ptlow"] = ptlow
result_dict["pthigh"] = pthigh
result_dict["etalow"] = etalow
result_dict["etahigh"] = etahigh


with open(outfile, 'wb') as f_dump:
    pickle.dump(result_dict, f_dump)