import awkward as ak
import numpy as np
import os
import glob
import ROOT


path = "/eos/user/l/ldellape/VBS/parquet"
category="baseline"
norm_luminosity=False
norm_events=False


#category="signal_AK8"
#category="signal_AK4"
#category="CR_TTbar"
#category="CR_L_Wjets"
#category="CR_R_Wjets"


if category == "baseline":
    parquet_patterns = [
        f"{path}/TTtoLNu2Q_HT-500_NJet-9_Hdamp-158_TuneCP5_13p6TeV_powheg-pythia8_2023_preBPix/{category}/*.parquet",
        f"{path}/ssWW_TT_mg5_madspin/{category}/*.parquet",
        f"{path}/ssWW_LL_mg5_madspin/{category}/*.parquet",
        f"{path}/WtoLNu-2Jets_PTLNu-100to200_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2023_preBPix/{category}/*.parquet"
    ]
    labels = ["TTtoLNu2Q", "ssWW_TT", "ssWW_LL", "WLtoLNu_2Jets"]
elif category=="signal_AK8":
    parquet_patterns = [
        f"{path}/ssWW_TT_mg5_madspin/SingleLepton_AK8/*.parquet",
        f"{path}/ssWW_LL_mg5_madspin/SingleLepton_AK8/*.parquet",
        f"{path}/TTtoLNu2Q_HT-500_NJet-9_Hdamp-158_TuneCP5_13p6TeV_powheg-pythia8_2023_preBPix/SingleLepton_AK8/*.parquet",
        f"{path}/WtoLNu-2Jets_PTLNu-100to200_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2023_preBPix/SingleLepton_AK8/*.parquet",
    ]
    labels = ["ssWW_TT", "ssWW_LL", "TTtoLNu2Q", "WLtoLNu_2Jets"]
elif category=="signal_AK4":
    parquet_patterns = [
        f"{path}/ssWW_TT_mg5_madspin/SingleLepton_AK4/*.parquet",
        f"{path}/ssWW_LL_mg5_madspin/SingleLepton_AK4/*.parquet",
        f"{path}/TTtoLNu2Q_HT-500_NJet-9_Hdamp-158_TuneCP5_13p6TeV_powheg-pythia8_2023_preBPix/SingleLepton_AK4/*.parquet",
        f"{path}/WtoLNu-2Jets_PTLNu-100to200_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2023_preBPix/SingleLepton_AK4/*.parquet",
    ]
    labels = ["ssWW_TT", "ssWW_LL", "TTtoLNu2Q", "WLtoLNu_2Jets"]
elif category=="CR_TTbar":
    parquet_patterns = [
        f"{path}/ssWW_TT_mg5_madspin/SingleEle_AK8_bjets_ttbar/*.parquet",
        f"{path}/ssWW_TT_mg5_madspin/SingleMuon_AK8_bjets_ttbar/*.parquet",
        f"{path}/ssWW_LL_mg5_madspin/SingleEle_AK8_bjets_ttbar/*.parquet",
        f"{path}/ssWW_LL_mg5_madspin/SingleMuon_AK8_bjets_ttbar/*.parquet",
        f"{path}/ssWW_TT_mg5_madspin/SingleEle_AK4_bjets_ttbar/*.parquet",
        f"{path}/ssWW_TT_mg5_madspin/SingleMuon_AK4_bjets_ttbar/*.parquet",
        f"{path}/ssWW_LL_mg5_madspin/SingleEle_AK4_bjets_ttbar/*.parquet",
        f"{path}/ssWW_LL_mg5_madspin/SingleMuon_AK4_bjets_ttbar/*.parquet",
        f"{path}/TTtoLNu2Q_HT-500_NJet-9_Hdamp-158_TuneCP5_13p6TeV_powheg-pythia8_2023_preBPix/SingleEle_AK8_bjets_ttbar/*.parquet",
        f"{path}/TTtoLNu2Q_HT-500_NJet-9_Hdamp-158_TuneCP5_13p6TeV_powheg-pythia8_2023_preBPix/SingleMuon_AK8_bjets_ttbar/*.parquet",
        f"{path}/WtoLNu-2Jets_PTLNu-100to200_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2023_preBPix/SingleEle_AK8_bjets_ttbar/*.parquet",
        f"{path}/WtoLNu-2Jets_PTLNu-100to200_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2023_preBPix/SingleMuon_AK8_bjets_ttbar/*.parquet"
    ]
    labels = ["ssWW_TT", "ssWW_TT", "ssWW_TT", "ssWW_TT", "ssWW_LL", "ssWW_LL", "ssWW_LL", "ssWW_LL","TTtoLNu2Q", "TTtoLNu2Q", "TTtoLNu2Q" , "TTtoLNu2Q", "WLtoLNu_2Jets","WLtoLNu_2Jets", "WLtoLNu_2Jets", "WtoLNu_2Jets"]
    

def load_awkward_parquet(pattern):
    files = glob.glob(pattern)
    arrays = [ak.from_parquet(f) for f in files]
    if len(arrays) == 0:
        print(f"Warning: no files found for pattern {pattern}")
        return None
    return ak.concatenate(arrays)


def sum_lumi(tag):
    if tag == "ssWW_LL":
        log_path = "log_ssWW_LL_mg5_madspin.txt"
    elif tag=="ssWW_TT":
        log_path = "log_ssWW_TT_mg5_madspin.txt"
    elif tag == "TTtoLNu2Q":
        log_path="log_TTtoLNu2Q.txt"
    elif tag == "WLtoLNu_2Jets":
        log_path ="log_WtoLNu-XJets.txt"
    else: 
        return 0.0
    total_lumi = 0.0

    if not os.path.exists(log_path):
        print(f"No log file found for dataset {dataset}")
        return 0.0

    with open(log_path, "r") as f:
        for line in f:
            parts = line.strip().split(',')
            for part in parts:
                if part.strip().startswith("lumi="):
                    lumi_value = float(part.split('=')[1])
                    total_lumi += lumi_value
    print(total_lumi)
    return total_lumi*1e-03

datasets = [load_awkward_parquet(p) for p in parquet_patterns]



lumis = [sum_lumi(label) for label in labels]
target_lumi = sum_lumi("ssWW_TT")
lumis = np.array(lumis, dtype=float)
print(lumis)
scaling_factors = target_lumi / lumis
for label, lumi, scale in zip(labels, lumis, scaling_factors):
    print(f"{label:<10} Lumi: {lumi:.3f}  Scale factor: {scale:.6f}")

c = ROOT.TCanvas("c", "Comparison", 800, 600)
ROOT.gStyle.SetOptStat(0)

field_name = ["lumi", "CleanFatJet_pt", "CleanFatJet_eta", "CleanFatJet_msoftdrop", "CleanFatJet_mass", "events_zepp_lep", "CleanSubJet_pair_zg", 
              "CleanFatJet_tau21", "events_nCleanJets", "events_nFatJet", "events_nCleanFatJets", "events_nBJetGood", "events_MT_lep_miss", "ElectronGood_pt", "MuonGood_pt", "V_dijet_candidate_mass", "V_dijet_candidate_pt",
              "VBS_dijet_system_pt", "VBS_dijet_system_mass", "VBS_dijet_system_deltaEta", "events_MT_ele_miss", "events_MT_mu_miss", "events_MT_lep_miss"]

number_events = []
for (data,label) in zip(datasets, labels):
    number_events.append(len(data))  
    print(f"{len(data)} events in {label}")  


for var in field_name:
    hists = []
    for i, (data, label) in enumerate(zip(datasets, labels)):
        if var in data.fields:
            vals = ak.flatten(data[var], axis=None)
            if var =="events_nCleanFatJets":
                print("ok")
                print(vals)
            if len(vals) == 0:
                print(f"[WARNING] No entries for variable '{var}' in dataset '{label}'")
                continue
            existing = next((h for (lab, h) in hists if lab == label), None)
            if existing:
                for v in vals:
                    existing.Fill(v)
            else:
                if var not in ["events_nCleanJet","events_nCleanFatJets", "events_nFatJet", "events_nJet", "events_nBJetGood"]:
                    h = ROOT.TH1F(f"h_{var}_{i}", label, 50, min(vals), max(vals))
                else:
                    h = ROOT.TH1F(f"h_{var}_{i}", label, 8, 0, 8)
                for v in vals:
                    h.Fill(v)
                h.SetLineColor(i + 1)
                h.SetLineWidth(2)
                if norm_luminosity is True:
                    h.Scale(scaling_factors[i])
                    
                hists.append((label, h))

    c = ROOT.TCanvas(f"c_{var}", f"Canvas for {var}", 800, 600)
    legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)

    for idx, (label, h) in enumerate(hists):
        draw_option = "hist" if idx == 0 else "hist same"
        h.Draw(draw_option)
        h.GetXaxis().SetTitle(var)
        h.GetYaxis().SetTitle("Entries")
        h.SetTitle(f"{var}   {category}")
        legend.AddEntry(h, label, "l")

    if hists:  
        hists[0][1].GetXaxis().SetTitle(var)
        hists[0][1].GetYaxis().SetTitle("Entries")

    legend.Draw()
    c.Update()
    c.Draw()
    os.makedirs(f"./plots/{category}", exist_ok=True)
    c.SaveAs(f"./plots/{category}/{var}_{category}.pdf")
