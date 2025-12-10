# Skim NANOAOD ntuples for L1 JEC/JERs derived using the MuonJet method
import ROOT
import time
import argparse

# Argument parser
parser = argparse.ArgumentParser(description='Preprocess L1 Scouting ntuples')
parser.add_argument('--input', '-i', type=str, help='Input file path')
parser.add_argument('--output', '-o', type=str, help='Output file path')
args = parser.parse_args()

fname = args.input
oname = args.output

print("Input file: {}".format(fname))
print("Output file: {}".format(oname))

# Start processing
start = time.time()

# Preprocess function call
ROOT.gInterpreter.Declare("""
    using Vbool = ROOT::RVec<bool>;
    using Vint = ROOT::RVec<int>;
    using Vfloat = ROOT::RVec<float>;

    auto getJetEt(Vfloat Jet_pt, Vfloat Jet_eta, Vfloat Jet_phi, Vfloat Jet_mass){
        Vfloat Jet_et(Jet_pt.size());
        for(size_t i = 0; i < Jet_pt.size(); ++i){
            ROOT::Math::PtEtaPhiMVector jet(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);
            Jet_et[i] = jet.Et();
        }
        return Jet_et;
    }

    // Get the matched reco jet eT
    auto getMatchedPt(Vfloat RecoJet_et, Vfloat RecoJet_eta, Vfloat RecoJet_phi, Vfloat Jet_pt, Vfloat Jet_eta, Vfloat Jet_phi) {
        Vfloat matchedPt(Jet_pt.size(), -1);
        Vint matchedRecoJetIdx(Jet_pt.size(), -1);

        for(size_t i = 0; i < Jet_pt.size(); ++i){
            float minDR = 9999;
            // Doing minDPT instead
            float minDPT = 9999;
            int bestRecoJetIdx = -1;

            for(size_t j = 0; j < RecoJet_et.size(); ++j){
                if (std::find(matchedRecoJetIdx.begin(), matchedRecoJetIdx.end(), j) != matchedRecoJetIdx.end()) continue;

                float dEta = Jet_eta[i] - RecoJet_eta[j];
                float dPhi = TVector2::Phi_mpi_pi(Jet_phi[i] - RecoJet_phi[j]);
                float dr = std::sqrt(dEta*dEta + dPhi*dPhi);
                float dpt = std::abs(Jet_pt[i] - RecoJet_et[j])/Jet_pt[i];
                
                //if((dr < minDR) && (dr < 0.2)){
                //    minDR = dr;
                //    matchedPt[i] = RecoJet_et[j];
                //    bestRecoJetIdx = j;
                //}

                // Using dPT matching instead of dR matching
                if((dpt < minDPT) && (dr < 0.2)){
                    minDPT = dpt;
                    matchedPt[i] = RecoJet_et[j];
                    bestRecoJetIdx = j;
                }
            }
            matchedRecoJetIdx[i] = bestRecoJetIdx;
        }
        return matchedPt;
    }
                          
    // Get the matched reco jet eta
    auto getMatchedEta(Vfloat RecoJet_et, Vfloat RecoJet_eta, Vfloat RecoJet_phi, Vfloat Jet_pt, Vfloat Jet_eta, Vfloat Jet_phi) {
        Vfloat matchedEta(Jet_pt.size(), -1);
        Vint matchedRecoJetIdx(Jet_pt.size(), -1);

        for(size_t i = 0; i < Jet_pt.size(); ++i){
            float minDR = 9999;
            // Doing minDPT instead
            float minDPT = 9999;
            int bestRecoJetIdx = -1;

            for(size_t j = 0; j < RecoJet_eta.size(); ++j){
                if (std::find(matchedRecoJetIdx.begin(), matchedRecoJetIdx.end(), j) != matchedRecoJetIdx.end()) continue;

                float dEta = Jet_eta[i] - RecoJet_eta[j];
                float dPhi = TVector2::Phi_mpi_pi(Jet_phi[i] - RecoJet_phi[j]);
                float dr = std::sqrt(dEta*dEta + dPhi*dPhi);
                float dpt = std::abs(Jet_pt[i] - RecoJet_et[j])/Jet_pt[i];

                //if((dr < minDR) && (dr < 0.2)){
                //    minDR = dr;
                //    matchedEta[i] = RecoJet_eta[j];
                //    bestRecoJetIdx = j;
                //}

                // Using dPT matching instead of dR matching
                if((dpt < minDPT) && (dr < 0.2)){
                    minDPT = dpt;
                    matchedEta[i] = RecoJet_eta[j];
                    bestRecoJetIdx = j;
                }
            }
            matchedRecoJetIdx[i] = bestRecoJetIdx;
        }
        return matchedEta;
    }

    // Get dR of nearest offline muon
    auto getdRNearestMuon(Vfloat Jet_eta, Vfloat Jet_phi, Vfloat Muon_eta, Vfloat Muon_phi) {
        Vfloat dRNearestMuon(Jet_eta.size(), 9999);
        for(size_t i = 0; i < Jet_eta.size(); ++i){
            for(size_t j = 0; j < Muon_eta.size(); ++j){
                float dEta = Jet_eta[i] - Muon_eta[j];
                float dPhi = TVector2::Phi_mpi_pi(Jet_phi[i] - Muon_phi[j]);
                float dr = std::sqrt(dEta*dEta + dPhi*dPhi);
                if (dr < dRNearestMuon[i]){
                    dRNearestMuon[i] = dr;
                }
            }
        }    
        return dRNearestMuon;                     
    }
""")

# Input data
data = ROOT.TChain()
data.Add("{}/Events".format(fname))
df = ROOT.RDataFrame(data)

# Event selection (skimming)
df = df.Filter("HLT_IsoMu24").Filter("nL1Jet > 1").Filter("nJet > 0")
# Muon pre-selections (akin to MuonJet JEC derivation)
df = df.Define('Muon_PassTightId','Muon_pfIsoId>=3&&Muon_mediumPromptId')
df = df.Define('goodmuonPt25','Muon_pt>25&&abs(Muon_pdgId)==13&&Muon_PassTightId')
df = df.Filter('Sum(goodmuonPt25)>=1','>=1 muon with p_{T}>25 GeV')
df = df.Define('badmuonPt10','Muon_pt>10&&abs(Muon_pdgId)==13&&Muon_PassTightId==0')
df = df.Filter('Sum(badmuonPt10)==0','No bad quality muon')

# RecoJet cleaning (optional)
#df = df.Define("absJetEta", "abs(Jet_eta)")
#df = df.Define("passPFJetID",
#    """
#    (absJetEta <= 2.6 && Jet_neHEF < 0.90 && Jet_neEmEF < 0.90 && Jet_nConstituents > 1 &&
#    Jet_muEF < 0.80 && Jet_chHEF > 0.01 && Jet_chMultiplicity > 0 && Jet_chEmEF < 0.80) ||
#
#    (absJetEta > 2.6 && absJetEta <= 2.7 && Jet_neHEF < 0.90 && Jet_neEmEF < 0.99 &&
#    Jet_muEF < 0.80 && Jet_chEmEF < 0.80) ||
#
#    (absJetEta > 2.7 && absJetEta <= 3.0 && Jet_neHEF < 0.9999) ||
#
#    (absJetEta > 3.0 && absJetEta <= 5.0 && Jet_neEmEF < 0.90 && Jet_neMultiplicity > 2)
#    """
#)
#df = df.Define('isCleanJet', 'passPFJetID && Jet_muEF < 0.5 && Jet_chEmEF < 0.5')
#df = df.Redefine('Jet_pt', 'Jet_pt[isCleanJet]')
#df = df.Redefine('Jet_eta', 'Jet_eta[isCleanJet]')
#df = df.Redefine('Jet_phi', 'Jet_phi[isCleanJet]')

# Sort jets by pT
df = df.Define( "L1Jet_ptorder", "Reverse(Argsort(L1Jet_pt))")
df = df.Redefine("L1Jet_pt", "Take(L1Jet_pt, L1Jet_ptorder)")
df = df.Redefine("L1Jet_eta", "Take(L1Jet_eta, L1Jet_ptorder)")
df = df.Redefine("L1Jet_phi", "Take(L1Jet_phi, L1Jet_ptorder)")

# Basic event selection
df = df.Filter("(L1Jet_pt[0] >= 20) && (L1Jet_pt[1] >= 20)")
df = df.Filter("std::abs(L1Jet_eta[0]) < 3.0 && std::abs(L1Jet_eta[1]) < 3.0")
df = df.Filter("Sum(L1Jet_pt == 1023.5) == 0", "Saturated L1 jet veto")

# Skim out all the jets
df = df.Define("L1Jet_skim", "(L1Jet_pt >= 20) && (abs(L1Jet_eta) < 3.0)")
df = df.Redefine("L1Jet_pt", "L1Jet_pt[L1Jet_skim]").Redefine("L1Jet_eta", "L1Jet_eta[L1Jet_skim]").Redefine("L1Jet_phi", "L1Jet_phi[L1Jet_skim]")

# Make sure that they're not dR matched to a muon
df = df.Define("L1Jet_dRNearestMuon", "getdRNearestMuon(L1Jet_eta, L1Jet_phi, Muon_eta, Muon_phi)")
df = df.Define("L1Jet_muonVeto", "L1Jet_dRNearestMuon > 0.4")
df = df.Redefine("L1Jet_pt", "L1Jet_pt[L1Jet_muonVeto]").Redefine("L1Jet_eta", "L1Jet_eta[L1Jet_muonVeto]").Redefine("L1Jet_phi", "L1Jet_phi[L1Jet_muonVeto]").Redefine("L1Jet_dRNearestMuon", "L1Jet_dRNearestMuon[L1Jet_muonVeto]")

# Define quantities
df = df.Define("Jet_et", "getJetEt(Jet_pt, Jet_eta, Jet_phi, Jet_mass)") 
df = df.Define("L1Jet_matchedRecoPt", "getMatchedPt(Jet_et, Jet_eta, Jet_phi, L1Jet_pt, L1Jet_eta, L1Jet_phi)")
df = df.Define("L1Jet_matchedRecoEta", "getMatchedEta(Jet_et, Jet_eta, Jet_phi, L1Jet_pt, L1Jet_eta, L1Jet_phi)")

# Filter L1 jets to only those that have matches
df = df.Define("L1Jet_isMatched", "L1Jet_matchedRecoPt > 0")
df = df.Define("MatchedJet_pt", "L1Jet_pt[L1Jet_isMatched]")
df = df.Define("MatchedJet_eta", "L1Jet_eta[L1Jet_isMatched]")
df = df.Define("MatchedJet_phi", "L1Jet_phi[L1Jet_isMatched]")
df = df.Define("MatchedJet_matchedRecoPt", "L1Jet_matchedRecoPt[L1Jet_isMatched]")
df = df.Define("MatchedJet_matchedRecoEta", "L1Jet_matchedRecoEta[L1Jet_isMatched]")
df = df.Define("MatchedJet_ptDiff", "MatchedJet_pt - MatchedJet_matchedRecoPt")
df = df.Define("MatchedJet_ptScale", "MatchedJet_pt / MatchedJet_matchedRecoPt")

# Save a snapshot with these quantities 
df.Snapshot(
    "Events", 
    "{}".format(oname), 
    ["MatchedJet_pt", "MatchedJet_eta", "MatchedJet_phi", "MatchedJet_matchedRecoPt", "MatchedJet_matchedRecoEta", "MatchedJet_ptDiff", "MatchedJet_ptScale"]    
)

# End processing
end = time.time()
print("Processing time: {:.2f} seconds".format(end - start))