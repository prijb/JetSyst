# Derive both scale and resolution from response plot
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

# Gaussian fit function
def gauss(x, a, mean, sig):
    return a * np.exp(-(x - mean)**2 / (2 * sig**2))

f = uproot.open(infile)
h_pt_response = f[f"pt_{etalow}to{etahigh}_scale_{etalow}to{etahigh}"].to_hist()

# Decode
etalow = decode_eta(etalow)
etahigh = decode_eta(etahigh)

h_pt_scale = h_pt_response[hist.loc(ptlow):hist.loc(pthigh), ::sum]
h_scale = h_pt_response[hist.loc(ptlow):hist.loc(pthigh):sum, :]

rebin_factor = 1
if rebin is not None: 
    #h_diff = h_diff[::(rebin*1j)]
    h_scale = h_scale[::(rebin*1j)]
    rebin_factor = rebin

# Get the mean pT for the scale and resolution calculation
pt_scale = h_pt_scale.axes[0].centers
pt_scale_values = h_pt_scale.values()
mean_pt_scale = np.average(pt_scale, weights=pt_scale_values)

x_scale = h_scale.axes[0].centers
y_scale = h_scale.values()
dy_scale = np.sqrt(y_scale)

####### Scale and resolution fit #########
center_pos_scale = np.argmax(y_scale)
x_scale_pos_low = np.max([center_pos_scale - (6 // rebin_factor), 0])
x_scale_pos_high = np.min([center_pos_scale + (6 // rebin_factor), len(x_scale)])
fit_x_scale = x_scale[x_scale_pos_low:x_scale_pos_high]
fit_y_scale = y_scale[x_scale_pos_low:x_scale_pos_high]
fit_dy_scale = dy_scale[x_scale_pos_low:x_scale_pos_high]
# Guess parameters
a_scale_guess = y_scale[center_pos_scale]
mean_scale_guess = np.average(fit_x_scale, weights=fit_y_scale)
sig_scale_guess = np.average((fit_x_scale - mean_scale_guess)**2, weights=fit_y_scale)**0.5

# Fit to Gaussian
try:
    p0_scale = [a_scale_guess, mean_scale_guess, sig_scale_guess]
    popt_scale, pcov_scale = curve_fit(gauss, fit_x_scale, fit_y_scale, p0=p0_scale, sigma=fit_dy_scale, absolute_sigma=True)
    #print(popt_scale)
    #print(pcov_scale)

    # Plot the fit
    plot_x_scale_gauss = np.linspace(fit_x_scale[0], fit_x_scale[-1], 100)
    plot_y_scale_gauss = gauss(plot_x_scale_gauss, *popt_scale)
    fig, ax = plt.subplots()
    hep.histplot(h_scale, ax=ax, label="Data", color="black", histtype="step", flow=None)
    ax.plot(plot_x_scale_gauss, plot_y_scale_gauss, label="Gaussian Fit", color="red")
    ax.set_xlabel("L1/RecoJet pT")
    ax.set_ylabel("Events")
    hep.cms.label(ax=ax, data=False, rlabel="Level-1 Scouting")
    ax.legend()
    ax.text(0.65, 0.75, f"pT bin [{ptlow}, {pthigh}]", transform=ax.transAxes, fontsize=16*1.2)
    ax.text(0.65, 0.70, f"eta bin [{etalow}, {etahigh}]", transform=ax.transAxes, fontsize=16*1.2)
    plt.savefig(f"{args.plotdir}/fit_scale_pt_{int(ptlow)}_{int(pthigh)}_eta_{encode_eta(etalow)}_{encode_eta(etahigh)}.png")

    ## Histogram stats
    mean_scale = np.average(fit_x_scale, weights=fit_y_scale)
    rms_scale = np.average((fit_x_scale - mean_scale)**2, weights=fit_y_scale)**0.5

    # Bootstrap to get the errors on the mean and RMS with fitting
    mean_scale_bootstrap = []
    rms_scale_bootstrap = []
    mean_scale_fit_bootstrap = []
    std_scale_fit_bootstrap = []

    for i_bootstrap in range(100):
        y_scale_bootstrap = np.random.poisson(y_scale)

        x_scale_pos_bootstrap = np.argmax(y_scale_bootstrap)
        x_scale_pos_low_bootstrap = np.max([x_scale_pos_bootstrap - (6 // rebin_factor), 0])
        x_scale_pos_high_bootstrap = np.min([x_scale_pos_bootstrap + (6 // rebin_factor), len(x_scale)])
        fit_x_scale_bootstrap = x_scale[x_scale_pos_low_bootstrap:x_scale_pos_high_bootstrap]
        fit_y_scale_bootstrap = y_scale_bootstrap[x_scale_pos_low_bootstrap:x_scale_pos_high_bootstrap]
        fit_dy_scale_bootstrap = dy_scale[x_scale_pos_low_bootstrap:x_scale_pos_high_bootstrap]

        a_scale_guess_bootstrap = y_scale_bootstrap[x_scale_pos_bootstrap]
        mean_scale_guess_bootstrap = np.average(fit_x_scale_bootstrap, weights=fit_y_scale_bootstrap)
        sig_scale_guess_bootstrap = np.average((fit_x_scale_bootstrap - mean_scale_guess_bootstrap)**2, weights=fit_y_scale_bootstrap)**0.5
        try:
            p0_scale_bootstrap = [a_scale_guess_bootstrap, mean_scale_guess_bootstrap, sig_scale_guess_bootstrap]
            popt_scale_bootstrap, pcov_scale_bootstrap = curve_fit(gauss, fit_x_scale_bootstrap, fit_y_scale_bootstrap, p0=p0_scale_bootstrap, sigma=fit_dy_scale_bootstrap, absolute_sigma=True)

            mean_scale_bootstrap.append(mean_scale_guess_bootstrap)
            rms_scale_bootstrap.append(sig_scale_guess_bootstrap)

            mean_scale_fit_bootstrap.append(popt_scale_bootstrap[1])
            std_scale_fit_bootstrap.append(popt_scale[2])

        except:
            print(f"Bootstrap fit {i_bootstrap} failed, skipping")

    mean_scale_bootstrap = np.array(mean_scale_bootstrap)
    rms_scale_bootstrap  = np.array(rms_scale_bootstrap)
    print(f"Number of successful bootstraps: {len(mean_scale_bootstrap)}/100")

    mean_scale_fit = popt_scale[1]
    std_scale_fit = popt_scale[2]
    mean_scale_fit_err = np.sqrt(pcov_scale[1, 1])
    std_scale_fit_err = np.sqrt(pcov_scale[2, 2])

    # From bootstrapping
    mean_scale_err = np.std(mean_scale_bootstrap)
    rms_scale_err = np.std(rms_scale_bootstrap)
except:
    # Plot the histogram without fit
    fig, ax = plt.subplots()
    hep.histplot(h_scale, ax=ax, label="Data", color="black", histtype="step", flow=None)
    ax.set_xlabel("L1/RecoJet pT")
    ax.set_ylabel("Events")
    hep.cms.label(ax=ax, data=False, rlabel="Level-1 Scouting")
    ax.legend()
    ax.text(0.65, 0.75, f"pT bin [{ptlow}, {pthigh}]", transform=ax.transAxes, fontsize=16*1.2)
    ax.text(0.65, 0.70, f"eta bin [{etalow}, {etahigh}]", transform=ax.transAxes, fontsize=16*1.2)
    plt.savefig(f"{args.plotdir}/fit_scale_pt_{int(ptlow)}_{int(pthigh)}_eta_{encode_eta(etalow)}_{encode_eta(etahigh)}.png")

    print(f"Scale fit failed for pt bin [{ptlow}, {pthigh}] and eta bin [{etalow}, {etahigh}]")

    ## Histogram stats
    mean_scale = np.average(fit_x_scale, weights=fit_y_scale)
    rms_scale = np.average((fit_x_scale - mean_scale)**2, weights=fit_y_scale)**0.5

 
    # Bootstrap to get the errors on the mean and RMS without fitting
    mean_scale_bootstrap = []
    rms_scale_bootstrap = []
    for i_bootstrap in range(100):
        y_scale_bootstrap = np.random.poisson(y_scale)

        x_scale_pos_bootstrap = np.argmax(y_scale_bootstrap)
        x_scale_pos_low_bootstrap = np.max([x_scale_pos_bootstrap - (6 // rebin_factor), 0])
        x_scale_pos_high_bootstrap = np.min([x_scale_pos_bootstrap + (6 // rebin_factor), len(x_scale)])
        fit_x_scale_bootstrap = x_scale[x_scale_pos_low_bootstrap:x_scale_pos_high_bootstrap]
        fit_y_scale_bootstrap = y_scale_bootstrap[x_scale_pos_low_bootstrap:x_scale_pos_high_bootstrap]
        fit_dy_scale_bootstrap = dy_scale[x_scale_pos_low_bootstrap:x_scale_pos_high_bootstrap]

        a_scale_guess_bootstrap = y_scale_bootstrap[x_scale_pos_bootstrap]
        mean_scale_guess_bootstrap = np.average(fit_x_scale_bootstrap, weights=fit_y_scale_bootstrap)
        sig_scale_guess_bootstrap = np.average((fit_x_scale_bootstrap - mean_scale_guess_bootstrap)**2, weights=fit_y_scale_bootstrap)**0.5

        mean_scale_bootstrap.append(mean_scale_guess_bootstrap)
        rms_scale_bootstrap.append(sig_scale_guess_bootstrap)

    mean_scale_bootstrap = np.array(mean_scale_bootstrap)
    rms_scale_bootstrap = np.array(rms_scale_bootstrap)
    print(f"Number of successful bootstraps: {len(mean_scale_bootstrap)}/100")

    mean_scale_fit = mean_scale
    std_scale_fit = rms_scale
    mean_scale_fit_err = np.std(mean_scale_bootstrap)
    std_scale_fit_err = np.std(rms_scale_bootstrap)

    mean_scale_err = np.std(mean_scale_bootstrap)
    rms_scale_err = np.std(rms_scale_bootstrap)

print("\nFit results")
print("Scale")
print(f"  Mean x: {mean_pt_scale:.2f}")
print(f"  Guess Mean: {mean_scale:.2f}")
print(f"  Fit Mean: {mean_scale_fit:.2f} +/- {mean_scale_fit_err:.2e}")
print("Resolution")
print(f"  Guess RMS: {rms_scale:.2f}")
print(f"  Fit Std: {std_scale_fit:.2f} +/- {std_scale_fit_err:.2e}")

import pickle
result_dict = {}
fit_dict = {}

fit_dict["x"] = mean_pt_scale
fit_dict["mean_guess"] = mean_scale
fit_dict["mean_guess_err"] = mean_scale_err
fit_dict["mean_fit"] = mean_scale_fit
fit_dict["mean_fit_err"] = mean_scale_fit_err
fit_dict["std_guess"] = rms_scale
fit_dict["std_guess_err"] = rms_scale_err
fit_dict["std_fit"] = std_scale_fit
fit_dict["std_fit_err"] = std_scale_fit_err 

result_dict["fit"] = fit_dict
# Add metadata
result_dict["ptlow"] = ptlow
result_dict["pthigh"] = pthigh
result_dict["etalow"] = etalow
result_dict["etahigh"] = etahigh

with open(outfile, 'wb') as f_dump:
    pickle.dump(result_dict, f_dump)