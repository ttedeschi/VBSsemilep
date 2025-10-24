import awkward as ak
import math
from pocket_coffea.lib.cut_definition import Cut
from pocket_coffea.lib.cut_functions import get_nObj_min, get_HLTsel, get_nPVgood, goldenJson, eventFlags, get_nElectron, get_nMuon, get_nObj_eq
from pocket_coffea.lib.objects import get_dilepton




#############################################################
# VBS TOPOLOGY                                              #
#############################################################
def VBS_topology(events, params, year, sample, **kwargs):
    mask = (
        (events.VBS_dijet_system.mass > params["mass"] )
        & 
        (events.VBS_dijet_system.deltaEta > params["deltaEta"])
        &
        (events.VBS_dijet_system.pt1 > params["pt_leadingJet"])
        &
        (
            (
                (events.nCleanFatJets == 1) &
                (events.nCleanJets >= params["nJet_with_FatJet"])
            )
            |
            (
                (events.nCleanJets >= params["nJet"]) &
                (events.nCleanFatJets == 0)
            )
        )
    )
    return ak.where(ak.is_none(mask), False, mask)
VBS_jets_presel = Cut(
    name="VBS_jets_presel",
    params = {
        "mass" : 400,
        "deltaEta" : 2.5,
        "pt_leadingJet" : 50,
        "nJet" : 4,
        "nJet_with_FatJet" : 2,
    },
    function=VBS_topology,
)
#############################################################
#############################################################



#############################################################
# semileptonic, requirements on leptons and/or MET          #
#############################################################
def semileptonic(events, params, year, sample, **kwargs):
    single_electron = events.nElectronGood == 1
    single_muon = events.nMuonGood == 1
    double_electron = events.nElectronGood == 2
    double_muon = events.nMuonGood == 2
    if params["W"] is True:
        mask = (
                (   single_electron
                  & (ak.firsts(events.LeptonGood.pt) > params["pt_leading_electron"])
                  & (events.MET.pt > params["met_electron"])
                )
                | 
                (
                   single_muon
                   & (ak.firsts(events.LeptonGood.pt) > params["pt_leading_muon"])
                   & (events.MET.pt > params["met_muon"])
                )
        )
    elif params["Z"] is True:
        muon_pt = ak.pad_none(events.MuonGood.pt, 2)
        muon_charge = ak.pad_none(events.MuonGood.charge, 2)
        electron_pt = ak.pad_none(events.ElectronGood.pt, 2)
        electron_charge = ak.pad_none(events.ElectronGood.charge, 2)
        mu_mask = (
                    (muon_pt[:, 0] > params["pt_leading_muon"]) &
                    (muon_pt[:, 1] > params["pt_subleading_muon"]) &
                    (muon_charge[:,0] != muon_charge[:, 1])
        )
        mu_mask = ak.fill_none(mu_mask, False)
        el_mask = (
            (electron_pt[:, 0] > params["pt_leading_electron"]) &
            (electron_pt[:, 1] > params["pt_subleading_electron"]) &
            (electron_charge[:,0] != electron_charge[:,1])
        )
        el_mask = ak.fill_none(el_mask, False)
        mask = mu_mask | el_mask
    return ak.where(ak.is_none(mask), False, mask)
semileptonic_preselW = Cut(
    name="semileptonic_preselW", 
    params = {
        "W" : True,
        "pt_leading_electron" : 30,
        "pt_leading_muon" : 30,
        "eta_max_lep" : 2.5, 
        "nJet" : 4, 
        "nFatJets" : 1,
        "nJet_with_FatJet" : 2,
        "met_electron" : 30,
        "met_muon" : 30,
    },
    function = semileptonic,
) 
semileptonic_preselZ = Cut(
    name="semileptonic_preselZ", 
    params = {
        "W" : False,
        "Z" : True,
        "pt_leading_electron" : 20,
        "pt_subleading_electron" : 20,
        "pt_leading_muon" : 20,
        "pt_subleading_muon" : 20,
        "nJet" : 1, 
        "nFatJets" : 1,
        "nJet_with_FatJet" : 1,
        "met" : 20,
    },
    function = semileptonic,
) 
semileptonic_preselDY = Cut(
    name="semileptonic_preselDY",
    params = {
        "pt_leading_electron" : 20,
        "pt_subleading_electron" : 20,
        "pt_leading_muon" : 20,
        "pt_subleading_muon" : 20,
        "nJet" : 1, 
        "nFatJets" : 1,
        "nJet_with_FatJet" : 1,
        "met" : 20,        
    },
    function = semileptonic,
)
#############################################################
#############################################################




#############################################################
# V jet mass from AK8 or 2AK4 (pair wuth the mass closer to W/Z mass)
#############################################################
def Vjet_mass(events, params, year, sample, **kwargs):
    fj_mask = (
    (
        (events.CleanFatJet.msoftdrop > params["mass_min"]) &
        (events.CleanFatJet.msoftdrop < params["mass_max"]) & 
        (events.CleanFatJet.tau21 < params["tau21"])
    )
    |
    (
        (events.V_dijet_candidate.mass > params["mass_min"]) &
        (events.V_dijet_candidate.mass < params["mass_max"]) 
    )
    )
    mask = ak.any(fj_mask, axis=1)
    return ak.where(ak.is_none(mask), False, mask)
Vjet_massZ = Cut(
    name="Vjet_massZ",
    params={
        "VV": True,
        "mass_min": 70,
        "mass_max": 110,
        
    },
    function=Vjet_mass,
)
Vjet_massW = Cut(
    name="Vjet_massW",
    params={
        "VV": True,
        "mass_min": 65,
        "mass_max": 105,
        "nJet_min" : 4,
        "tau21":0.45,
    },
    function=Vjet_mass,
)
#############################################################
#############################################################



#############################################################
# number of b-tagged (central) jets (for ttbar etc.)        #
#############################################################
def Bjets_presel(events, params, year, sample, **kwargs):
    b_mask = (events["nBJetGood"] >= params["nBjets"])
    return ak.where(ak.is_none(b_mask), False, b_mask)
Bjets_presel = Cut(
    name="Bjets_presel",
    params={
        "nBjets" : 1,
    },
    function=Bjets_presel,
)
#############################################################
#############################################################



#################################################################
# v-jet mass outside Vjet_mass (contamination from W+jets etc.) #           
#################################################################
# this always return true but split the sample
def Vjet_massSide(events, params, year, sample, **kwargs):
    if params.get("side", True):
        left_band_jets = (
            (   (ak.any(events.CleanFatJet.msoftdrop < params["mass_min"], axis=1)) &
                (events.nCleanFatJets == 1)
            )
            |
            (
                (events.V_dijet_candidate.mass < params["mass_min"]) &
                (events.nCleanJets >= params["nJet_min"])
            )
        )
        right_band_jets = (
            (
                (ak.any(events.CleanFatJet.msoftdrop > params["mass_max"],axis=1)) &
                (events.nCleanFatJets == 1)
            )
            |
            (
                (events.V_dijet_candidate.mass > params["mass_max"]) &
                (events.nCleanJets >= params["nJet_min"])
            )
        )                
        mask = left_band_jets | right_band_jets
        if params["W"] is True:
            events["Wjets_sideL"] = left_band_jets
            events["Wjets_sideR"] = right_band_jets
        else: 
            events["Zjets_sideL"] = left_band_jets
            events["Zjets_sideR"] = right_band_jets
        return ak.ones_like(mask, dtype=bool)
Vjet_massSideW = Cut(
    name="Vjet_massSide",
    params={
        "W" : True,
        "mass_min" : 60,
        "mass_max" : 110,
        "nJet_min" : 4,
    },
    function=Vjet_massSide,
)
Vjet_massSideZ = Cut(
    name="Vjet_massSide",
    params={
        "W" :  False,
        "mass_min" : 70,
        "mass_max" : 120,
        "nJet_min" : 4,
    },
    function=Vjet_massSide,
)
def sideL(events, params, year, sample, **kwargs):
    if params["boson"] =="W":
        mask = (events.Wjets_sideL == True)
    elif params["boson"]=="Z":
        mask = (events.Zjets_sideL == True)
    return ak.where(ak.is_none(mask), False, mask)
def sideR(events, params, year, sample, **kwargs):
    if params["boson"] == "W":
        mask = (events.Wjets_sideR == True)
    elif params["boson"] =="Z":
        mask = (events.Zjets_sideL == True)
    return ak.where(ak.is_none(mask), False, mask)
Wjet_sideL = Cut(
    name="Wjet_sideL",
    params = {"boson": "W"},
    function=sideL,
)
Wjet_sideR = Cut(
    name="Wjet_sideR",
    params={"boson" : "W"},
    function=sideR,
)
Zjet_sideL = Cut(
    name="Zjet_sideL",
    params={"boson" : "Z"},
    function=sideL,
)
Zjet_sideR = Cut(
    name="Zjet_sideL",
    params={"boson": "Z"},
    function=sideR,
)
#############################################################
#############################################################




#############################################################
# transverse mass requirements form W->lv                   #
#############################################################
def Wtransverse_mass(events, params, year, sample, **kwargs):
    mask = (
        (events.MT_lep_miss < params["transverse_max"]) &
        (events.MT_lep_miss > 0)
    )
    return ak.where(ak.is_none(mask), False, mask)
Wtransverse_mass_presel = Cut(
    name="Wtransverse_mass_presel",
    params={
        "transverse_max" : 180,
    },
    function=Wtransverse_mass,
)
#############################################################
#############################################################





#############################################################
# Z->ll dilepton mass                                       #
#############################################################
def dilepton_mass(events, params, year, sample, **kwargs):
    double_muon = events.nMuonGood == 2
    double_electron = events.nElectronGood == 2
    if params["VV"] is True:
        mask = (
            ( 
                (double_muon) &
                (events.dimuon_candidate.mass > params["mass_min"]) & 
                (events.dimuon_candidate.mass < params["mass_max"])
            )
            |
            ( 
                (double_electron) &
                (events.dielectron_candidate.mass > params["mass_min"]) &
                (events.dielectron_candidate.mass < params["mass_max"])
            )
        )
    elif params["DY"] is True:
        mask = (
            (events.dilepton_candidate.mass > params["mass_max"] | events.dilepton_candidate.mass < params["mass_min"])
        )
    return ak.where(ak.is_none(mask), False, mask)

dilepton_massZ = Cut(
    name="dilepton_massZ",
    params = {
        "VV" : True,
        "mass_min" : 70,
        "mass_max" : 110,
    },
    function= dilepton_mass,
)
dilepton_massDY = Cut(
    name="dilepton_massDY",
    params = {
        "DY" : True,
        "mass_min" : 70,
        "mass_max" : 110,
    },
    function = dilepton_mass,
)   
#############################################################
#############################################################





#############################################################
# leptons flavour (Z->ll vs tt)                             #
#############################################################
def check_flavour(events, params, year, sample, **kwargs):
    double_leptons = events.nLeptonGood == 2
    single_electron = events.nElectronGood == 1
    single_muon = events.nMuonGood == 1
    if params["VV"] is True:
        mask=(
            double_leptons
            &
            (
                (single_electron & ~single_muon)
                |
                (single_muon & ~single_electron)
            )
        )
    else: 
        mask=(
            double_leptons
            &
            single_electron
            & 
            single_muon
        )
    return ak.where(ak.is_none(mask), False, mask)
check_flavour_SF = Cut(
    name="check_flavour_SF",
    params={ "VV" : True, },
    function=check_flavour,
)
check_flavour_OF =Cut(
    name="check_flavour_OF",
    params={ "VV":False,},
    function=check_flavour,
)
#############################################################
#############################################################













def cut_function(events, params, year, sample, **kwargs):
    masks = []
    for i, c in enumerate(params["coll"]):
        objs = getattr(events, c)
        ptmin = params["minpt"][i]
        nmin = params["nMin"][i]

        if ptmin is not None:
            objs = objs[objs.pt > ptmin]

        mask = ak.num(objs) >= nmin
        masks.append(mask)
    combined_mask = masks[0]
    for m in masks[1:]:
        combined_mask = combined_mask | m
    return ak.where(ak.is_none(combined_mask), False, combined_mask)

def cut_function_2(events, params, year, sample, **kwargs):
    masks = []
    for i, c in enumerate(params["coll"]):
        objs = getattr(events, c)
        nmin = params["nMin"][i]
        mask = ak.num(objs) == nmin
        masks.append(mask)
    combined_mask = masks[0]
    for m in masks[1:]:
        combined_mask = combined_mask | m

    return ak.where(ak.is_none(combined_mask), False, combined_mask)


def skim_single_leptons(events, params, year, sample, **kwargs):
    mask = (
        (ak.num(events.Muon) >=1)
        |
        (ak.num(events.Electron) >=1)
    )
    return ak.where(ak.is_none(mask), False, mask)
skim_sigle_lepton= Cut(
    name="skim_single_leptons",
    params={},
    function=skim_single_leptons
)

def skim_double_leptons(events, params, year, sample, **kwargs):
    mask = (
        (ak.num(events.Muon) >= 2)
        |
        (ak.num(events.Electron) >= 2)
    )
    return ak.where(ak.is_none(mask), False, mask)

skim_double_lepton = Cut(
    name="skim_double_leptons",
    params={},
    function=skim_double_leptons,
)




def get_nObj_min_or(nMin, minpt,coll, name=None):
    if not (len(coll) == len(nMin) and (minpt is None or len(minpt) == len(coll))):
        raise ValueError("coll, nMin, and minpt (if provided) must have same length")

    if name is None:
        name = "cut_" + "_or_".join([
            f"{c}_min{nMin[i]}" + (f"_pt{minpt[i]}" if minpt else "")
            for i, c in enumerate(coll)
        ])
    return Cut(
        name=name,
        params={
            "coll": coll,
            "nMin": nMin,
            "minpt": minpt if minpt else [None] * len(coll),
        },
        function=cut_function,
    )

def single_good_electron(events, params, year, sample, **kwargs):
    mask = (
        ak.num(events.ElectronGood) >= 2
        )
    return ak.where(ak.is_none(mask), False, mask)  

def single_good_muon(events, params, year, sample, **kwargs):
    mask = (
        ak.num(events.MuonGood) >= 2
    )  
    return ak.where(ak.is_none(mask), False, mask)
def single_good_lepton(events, params, year, sample, **kwargs):
    ele_mask = single_good_electron(events, params, year, sample, **kwargs)
    muon_mask = single_good_muon(events, params, year, sample, **kwargs)
    combined_mask = ele_mask & muon_mask
    return ak.where(ak.is_none(combined_mask), False, combined_mask)   












SingleEle = Cut(
    name="SingleGoodEle",
    params = {},
    function=single_good_electron,
)
SingleMuon = Cut(
    name="SingleGoodMuon", 
    params = {},
    function=single_good_muon,
)
SingleLepton = Cut(
    name="SingleGoodLeptons",
    params={},
    function=single_good_lepton,
)
def get_nObj_eq_or(nMin, coll, name=None):
    if name is None:
        name = "cut_" + "_or_".join([
            f"{c}_min{nMin[i]}"
            for i, c in enumerate(coll)
        ])
    return Cut(
        name=name,
        params={
            "coll": coll,
            "nMin": nMin,
        },
        function=cut_function_2,
    )