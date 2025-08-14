from custom_cut_functions import *
from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.cut_definition import Cut
from pocket_coffea.lib.cut_functions import get_nObj_min, get_HLTsel, get_nPVgood, goldenJson, eventFlags, get_nElectron, get_nMuon, get_nObj_eq, count_objects_eq
from pocket_coffea.parameters.cuts import passthrough
from pocket_coffea.parameters.histograms import *
import workflowVBS
from workflowVBS import VBS_WV_Processor
from skim_dictionary import skim_dict
from pocket_coffea.lib.weights.common import common_weights
import cloudpickle
from pocket_coffea.lib.columns_manager import ColOut

import custom_cut_functions 

cloudpickle.register_pickle_by_value(workflowVBS)
cloudpickle.register_pickle_by_value(custom_cut_functions)

import os
localdir = os.path.dirname(os.path.abspath(__file__))

from pocket_coffea.parameters import defaults
default_parameters = defaults.get_default_parameters()
defaults.register_configuration_dir("config_dir", localdir+"/parameters")


parameters = defaults.merge_parameters_from_files(default_parameters,
                                    f"{localdir}/parameters/object_presel.yaml",
                                    f"{localdir}/parameters/btagging.yaml",
                                    update=True                                    
                                    )



parquet = False
cfg = Configurator(
    parameters=parameters,
    datasets = {
        "tag" : "ZZ",
        "jsons" : [#f"{localdir}/datasets/VBS_mc_132X_Summer23wmLHEGS.json",
                   f"{localdir}/datasets/TTbar/TTtoLNu2Q_HT-500_NJet-9_Hdamp-158_TuneCP5_13p6TeV_powheg-pythia8.json",
                   f"{localdir}/datasets/WJets/WtoLNu-2Jets_PTLNu-100to200_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8.json",
                   #f"{localdir}/datasets/TbarWplus/*.json"
                   ],
        "filter" : {
            "samples" : ["ssWWTT", "ssWWLL", "TTtoLNu2Q_HT-500_NJet-9_Hdamp-158_TuneCP5_13p6TeV_powheg-pythia8", "WtoLNu-2Jets_PTLNu-100to200_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8"],
        }
        }, 
    workflow=VBS_WV_Processor,
    skim = [
            get_nPVgood(1), 
            get_nObj_min_or(
                [skim_dict["ZZ"]["Muon"]["ptMin"],skim_dict["ZZ"]["Electron"]["ptMin"]], 
                [skim_dict["ZZ"]["Muon"]["nMin"], skim_dict["ZZ"]["Electron"]["nMin"]],
                ["Muon", "Electron"]),
            ],
    
    # signal
    # Vjet_massSide always return true, flag for Wjets signal and bkg, to be moved in the categorization
    preselections=[VBS_jets_presel, semileptonic_preselZ, Vjet_massSideZ],

    
    
    #preselections = [passthrough],
    # AK8 --> V boson from FatJet
    # AK4 --> V boson from 2 Jets
    
    
    # Vjet_massZ ---> on-shell W bosons
    # Wjet_sideX ---> off-shell W bosons
    categories= {
        "baseline" : [passthrough],
        
        
        "DoubleEle_AK8" : [get_nElectron(2, coll="ElectronGood"), get_nObj_eq(1, coll="CleanFatJet"), get_nObj_eq(0, coll="BJetGood"), Vjet_massZ, dilepton_massZ],
        "DoubleEle_AK4" : [get_nElectron(2, coll="ElectronGood"), get_nObj_min(4 , coll="CleanJet"),  get_nObj_eq(0, coll="CleanFatJet"), get_nObj_eq(0, coll="BJetGood"), Vjet_massZ, dilepton_massZ],
        "DoubleMuon_AK8" : [get_nMuon(2, coll="MuonGood"), get_nObj_eq(1, coll="CleanFatJet"), get_nObj_eq(0, coll="BJetGood"), Vjet_massZ, dilepton_massZ],
        "DoubleMuon_AK4" : [get_nMuon(2, coll="MuonGood"), get_nObj_min(4 , coll="CleanJet"),  get_nObj_eq(0, coll="CleanFatJet"), get_nObj_eq(0, coll="BJetGood"), Vjet_massZ, dilepton_massZ],
        "DoubleLepton_AK8" : [get_nObj_eq(2, coll="LeptonGood"), get_nObj_eq(1, coll="CleanFatJet") , get_nObj_eq(0, coll="BJetGood"), Vjet_massZ, dilepton_massZ],
        "DoubleLepton_AK4" : [get_nObj_eq(2, coll="LeptonGood"), get_nObj_min(4, coll="CleanJet") , get_nObj_eq(0, coll="CleanFatJet"), get_nObj_eq(0, coll="BJetGood"), Vjet_massZ, dilepton_massZ],

        # ttbar/tbarWplus
        "DoubleEle_AK8_bjets_ttbar" : [get_nElectron(1, coll="ElectronGood"), get_nObj_eq(1 , coll="CleanFatJet"), get_nObj_min(2, coll="BJetGood")],
        "DoubleEle_AK4_bjets_ttbar" : [get_nElectron(1, coll="ElectronGood"), get_nObj_min(4 , coll="CleanJet"), get_nObj_min(2, coll="BJetGood")],
        "DoubleMuon_AK8_bjets_ttbar" : [get_nMuon(1, coll="MuonGood"), get_nObj_eq(1 , coll="CleanFatJet"), get_nObj_eq(2, coll="BJetGood")],
        "DoubleMuon_AK4_bjets_ttbar" : [get_nMuon(1, coll="MuonGood"), get_nObj_min(4, coll="CleanJet"), get_nObj_min(2, coll="BJetGood")],
        "DoubleLepton_AK8_bjets_ttbar" : [get_nObj_eq(1, coll="LeptonGood"), get_nObj_min(4, coll="CleanJet"), get_nObj_min(2, coll="BJetGood")],
        "DoubleLepton_AK4_bjets_ttbar" : [get_nObj_eq(1, coll="LeptonGood"), get_nObj_min(4, coll="CleanJet"), get_nObj_min(2, coll="BJetGood")],
        
        #WtoLNu-XJets (check contamination in the jet mass window of SR)
        "DoubleEle_AK8_AK4_sideL_Wjets" : [get_nElectron(1, coll="ElectronGood"), Zjet_sideL],
        "DoubleEle_AK8_AK4_sideR_Wjets" : [get_nElectron(1, coll="ElectronGood"), Zjet_sideR],
        "DoubleMuon_AK8_AK4_sideL_Wjets" : [get_nMuon(1, coll="MuonGood"), Zjet_sideL],
        "DoubleMuon_AK8_AK4_sideR_Wjets" : [get_nMuon(1, coll="MuonGood"), Zjet_sideR],
        "DoubleLepton_AK8_AK4_sideL_Wjets" : [get_nObj_eq(1, coll="LeptonGood"), Zjet_sideL],
        "DoubleLepton_AK8_AK4_sideR_Wjets" : [get_nObj_eq(1, coll="LeptonGood"), Zjet_sideR]
    },    
    weights_classes = common_weights,
    weights = {
        "common": {
            "inclusive": ["genWeight","lumi","XS",
                          "pileup",
                          ],
            "bycategory" : {
            }
        },
        "bysample": {
        }
    },
    variations = {
        "weights": {
            "common": {
                "inclusive": [ 
                              ],
                "bycategory" : {
                }
            },
        "bysample": {
        }    
        },
    },
    variables = {},
    workflow_options = {
        "dump_columns_as_arrays_per_chunk": "root://eosuser.cern.ch//eos/user/l/ldellape/VBS/parquet/"
    },
    columns = {
        "common" : {
            "inclusive" : [ColOut("ElectronGood", ["pt","eta", "phi"], flatten=False),
                           ColOut("MuonGood", ["pt", "eta", "phi"], flatten=False),
                           ColOut("LeptonGood", ["pt", "eta", "phi"], flatten=False),
                           ColOut("CleanFatJet", ["pt", "eta", "phi", "tau1", "tau2", "tau21", "msoftdrop", "mass"], flatten=False),
                           ColOut("CleanJet", ["pt", "eta", "phi"], flatten=False),
                           ColOut("CleanJet_noVBS", ["pt", "eta", "phi"], flatten=False),
                           ColOut("VBS_dijet_system", ["mass", "pt", "deltaEta"], flatten=False),    
                           ColOut("V_dijet_candidate", ["mass", "pt", "deltaEta"], flatten=False),  
                           ColOut("events", ["zepp_ele", "zepp_muon", "zepp_lep"], flatten=False),
                           ColOut("events", ["MT_mu_miss", "MT_ele_miss", "MT_lep_miss"], flatten=False),
                           ColOut("events", ["nCleanJets", "nCleanFatJets", "nBJetGood", "nJet", "nFatJet"], flatten=False),
                           ColOut("CleanSubJet_pair", ["zg"], flatten=False),
                           ],
            "bycategory" : {},
            },            
    },   
)