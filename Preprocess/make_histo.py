# Makes histograms of the pT scale and difference for the matched jets
import ROOT
import time
import argparse

# Argument parser
parser = argparse.ArgumentParser(description='Preprocess L1 Scouting ntuples')
parser.add_argument('--input', '-i', type=str, help='Input file path')
parser.add_argument('--output', '-o', type=str, help='Output file path')
parser.add_argument('--ptrange', nargs='+', type=float, help='pT range for histograms (min,max)')
parser.add_argument('--etarange', nargs='+', type=float, help='Number of eta bins for histograms')
args = parser.parse_args()

fname = args.input
oname = args.output
ptrange = args.ptrange
pt_min = ptrange[0]
pt_max = ptrange[1]
etarange = args.etarange
eta_min = etarange[0]
eta_max = etarange[1]

print("Input file: {}".format(fname))
print("Output file: {}".format(oname))
print(f"pT range: {pt_min} to {pt_max}")
print(f"Eta range: {eta_min} to {eta_max}")

# Start processing
start = time.time()
ROOT.EnableImplicitMT(8)

# Input data
data = ROOT.TChain()
data.Add("{}/Events".format(fname))
df = ROOT.RDataFrame(data)

# Define plot quantities
# Setting scale cuts in terms of the L1 quantities 
df = df.Define("pt_cut_scale", f"(MatchedJet_pt >= {pt_min}) && (MatchedJet_pt < {pt_max})")
df = df.Define("eta_cut_scale", f"(MatchedJet_eta >= {eta_min}) && (MatchedJet_eta < {eta_max})")
df = df.Define("pt_eta_cut_scale", "pt_cut_scale && eta_cut_scale")
# Resolution cuts in terms of L1 quantities
df = df.Define("pt_cut_reso", f"(MatchedJet_pt >= {pt_min}) && (MatchedJet_pt < {pt_max})")
df = df.Define("eta_cut_reso", f"(MatchedJet_eta >= {eta_min}) && (MatchedJet_eta < {eta_max})")
df = df.Define("pt_eta_cut_reso", "pt_cut_reso && eta_cut_reso")
df = df.Define("plot_pt_scale", "MatchedJet_pt[pt_eta_cut_scale]") 
df = df.Define("plot_pt_reso", "MatchedJet_pt[pt_eta_cut_reso]")
df = df.Define("plot_scale", "MatchedJet_ptScale[pt_eta_cut_scale]")
df = df.Define("plot_diff", "MatchedJet_ptDiff[pt_eta_cut_reso]")

# Create histograms
#h_pt_scale = df.Histo1D(("h_pt_scale", "L1 Jet pT; pT (GeV); Events", 120, 0, 600), "plot_pt_scale") # If x axis is reco pT
h_pt_scale = df.Histo1D(("h_pt_scale", "L1 Jet pT; pT (GeV); Events", 30, pt_min, pt_max), "plot_pt_scale") # If x axis is L1 pT
h_pt_reso = df.Histo1D(("h_pt_reso", "L1 Jet pT; pT (GeV); Events", 30, pt_min, pt_max), "plot_pt_reso")
h_diff = df.Histo1D(("h_diff", "L1 - RecoJet pT; pT Difference (GeV); Events", 100, -100, 100), "plot_diff")
h_scale = df.Histo1D(("h_scale", "L1/RecoJet pT; Scale; Events", 50, 0, 3.0), "plot_scale")

f_out = ROOT.TFile(oname, "RECREATE")
h_pt_scale.Write()
h_pt_reso.Write()
h_diff.Write()
h_scale.Write()
f_out.Close()

end = time.time()
print(f"Processing completed in {end - start:.2f} seconds")