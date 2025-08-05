import awkward as ak
import math
from pocket_coffea.lib.cut_definition import Cut
from pocket_coffea.lib.cut_functions import get_nObj_min, get_HLTsel, get_nPVgood, goldenJson, eventFlags, get_nElectron, get_nMuon, get_nObj_eq
from pocket_coffea.lib.objects import get_dilepton

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

def VBS_jets(events, params, year, sample, **kwargs):
    mask = (
        (events.VBS_dijet_system.mass > params["mass"] )
        & 
        (events.VBS_dijet_system.deltaEta > params["deltaEta"])
    )
    return ak.where(ak.is_none(mask), False, mask)

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





def Bjets_presel(events, params, year, sample, **kwargs):
    print(f"b jets : {events.nBJetGood}")
    
    b_mask = (
        (
            (events.nBJetGood >= params["nBjets"]) &
            (events.nCleanFatJets == params["nFatJet"]) &
            (events.nCleanJets >= params["nCleanJet_with_FatJet"])
        )
        |
        (
            (events.nBJetGood >= params["nBjets"]) &
            (events.nCleanJets >= params["nJets"])
        )
    )
    
    return ak.where(ak.is_none(b_mask), False, b_mask)

Bjets_presel = Cut(
    name="Bjets_presel",
    params={
        "nBjets" : 2,
        "nFatJet" : 1,
        "nCleanJet_with_FatJet" : 2,
        "nJets" : 4,
    },
    function=Bjets_presel,
)



# 1 lepton + 4 jets OR 1 lepton + 1 FatJet + 2 Jets
def semileptonic(events, params, year, sample, **kwargs):
    print(f" clean fatjet pt: {events.CleanFatJet.pt}")
    single_electron = events.nElectronGood == 1
    single_muon = events.nMuonGood == 1
    double_electron = events.nElectronGood == 2
    double_muon = events.nMuonGood == 2
    print(f" MET: {events.MET.pt}")
    print(f" Muon: {events.MuonGood.pt}")
    print(f" Electron: {events.ElectronGood.pt}")
    if params["W"] is True:
        print("ok W")
        mask = (
            (events.nLeptonGood == 1)
        &
        (
            (
                single_electron
                & (
                    ak.firsts(events.LeptonGood.pt)
                    > params["pt_leading_electron"]
                )
                & (
                    ak.firsts(events.LeptonGood.eta) < params["eta_max_lep"]
                )
                & (
                    events.MET.pt > params["met_electron"]
                )
            )
            | (
                single_muon
                & (ak.firsts(events.LeptonGood.pt) > params["pt_leading_muon"])
                & (ak.firsts(events.LeptonGood.eta) < params["eta_max_lep"])
                & (events.MET.pt > params["met_muon"])
            )
        )
        & (events.nCleanJets >= params["nJet"])
        | (
            (events.nCleanFatJets >= params["nFatJets"])
            & (events.nCleanJets >= params["nJet_with_FatJet"])
          )
        )

    elif params["Z"] is True:
        mask = (
            (events.nLeptonGood == 2)
            & (
                (
                    double_muon
                    & (events.LeptonGood.pt[:, 0] > params["pt_leading_muon"])
                    & (events.LeptonGood.pt[:, 1] > params["pt_subleading_muon"])
                )
                | (
                    double_electron
                    & (events.LeptonGood.pt[:, 0] > params["pt_leading_electron"])
                    & (events.LeptonGood.pt[:, 1] > params["pt_subleading_electron"])
                )
            )
            & (
                (events.nCleanJets >= params["nJet"])
                | (
                    (events.nCleanFatJets >= params["nFatJets"])
                    & (events.nCleanJets >= params["nJet_with_FatJet"])
                )
            )
        )

    print(mask)

    return ak.where(ak.is_none(mask), False, mask)



# ------------- dilepton mass ------------------- #
def dilepton_mass(events, params, year, sample, **kwargs):
    if events.nElectronGood==2 and params["VV"] is True:
        dilepton_candidate = get_dilepton(events["ElectronGood"])
        mask = (
            dilepton_candidate.pt > params["mass_min"] & dilepton_candidate.mass < params["mass_max"]
        )
    elif events.nMuonGood==2 and params["VV"] is True:
        dilepton_candidate = get_dilepton(events["MuonGood"])
        mask = (
            dilepton_candidate.mass > params["mass_min"]  & dilepton_candidate.mass < params["mass_max"]
        )
    elif events.nElectronGood==2 and params["DY"] is True:
        dilepton_candidate = get_dilepton(events["ElectronGood"])
        mask = (
            dilepton_candidate.mass > params["mass_max"] & dilepton_candidate.mass < params["mass_min"]
        )
    elif events.nMuonGood==2 and params["DY"] is True:
        dilepton_candidate = get_dilepton(events["MuonGood"])
        mask = (
            dilepton_candidate.mass > params["mass_max"] & dilepton_candidate.mass < params["mass_min"]
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
dilepton_massW = Cut(
    name = "dilepton_massW",
    params = {
        "VV" : True,
        "mass_min" : 70,
        "mass_max" : 110,
    },
    function = dilepton_mass,
) 
# --------------------------------------------------------------- #



# ---------- FatJet mass or dijet mass -------------------------- #
def Vjet_mass(events, params, year, sample, **kwargs):
    print(f" again fj: {events.CleanFatJet.mass}")
    print(f" mass : {events.CleanFatJet.mass}")
    fj_mask = (
    (
        (events.CleanFatJet.mass > params["mass_min"]) &
        (events.CleanFatJet.mass < params["mass_max"]) 
    )
    |
    (
           (events.V_dijet_candidate.mass > params["mass_min"]) &
           (events.V_dijet_candidate.mass < params["mass_max"]) &
           (events.nCleanJets >= params["nJet_min"])
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
        "mass_min": 60,
        "mass_max": 120,
        "nFatJet" : 1,
        "nJet_min" : 4,
    },
    function=Vjet_mass,
)

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
        "met_electron" : 90,
        "met_muon" : 40,
    },
    function = semileptonic,
) 
semileptonic_preselZ = Cut(
    name="semileptonic_preselZ", 
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
VBS_jets_presel = Cut(
    name="VBS_jets_presel",
    params = {
        "mass" : 400,
        "deltaEta" : 2.5,
    },
    function=VBS_jets,
)
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
