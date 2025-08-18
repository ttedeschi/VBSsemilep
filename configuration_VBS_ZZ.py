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
        "jsons" : [f"{localdir}/datasets/VBS_mc_132X_Summer23wmLHEGS.json",
                   #f"{localdir}/datasets/TTto2L2Nu-2Jets/TTto2L2Nu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_redirector.json",
                   #f"{localdir}/datasets/DY/DYto2L-2Jets_MLL-10to50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8.json",
                   #f"{localdir}/datasets/DY/DYto2L-2Jets_MLL-50_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8.json".
                   #f"{localdir}/datasets/DY/DYto2L-4Jets_MLL-50to120_HT-400to800_TuneCP5_13p6TeV_madgraphMLM-pythia8.json",
                   ],
        "filter" : {
            "samples" : ["ZZTT", "ZZLL", "DYto2L-2Jets_MLL-10to50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8", "DYto2L-4Jets_MLL-50to120_HT-400to800_TuneCP5_13p6TeV_madgraphMLM-pythia8", "DYto2L-2Jets_MLL-50_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8", "TTto2L2Nu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8"],
        }
        }, 
    workflow=VBS_WV_Processor,
    skim = [
            get_nPVgood(1), 
            skim_double_lepton,
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
        
        # signal
        "DoubleEle_AK8" : [get_nElectron(2, coll="ElectronGood"), get_nObj_eq(1, coll="CleanFatJet"), Vjet_massZ, dilepton_massZ],
        "DoubleEle_AK4" : [get_nElectron(2, coll="ElectronGood"), get_nObj_min(4 , coll="CleanJet"),  get_nObj_eq(0, coll="CleanFatJet"), Vjet_massZ, dilepton_massZ],
        "DoubleMuon_AK8" : [get_nMuon(2, coll="MuonGood"), get_nObj_eq(1, coll="CleanFatJet"), Vjet_massZ, dilepton_massZ],
        "DoubleMuon_AK4" : [get_nMuon(2, coll="MuonGood"), get_nObj_min(4 , coll="CleanJet"),  get_nObj_eq(0, coll="CleanFatJet"), Vjet_massZ, dilepton_massZ],

        # TTto2L2Nu-2Jets
        "AK8_OF" : [get_nMuon(1, coll="MuonGood"), get_nElectron(1, coll="ElectronGood"), get_nObj_eq(1, coll="CleanFatJet"), Vjet_massZ, dilepton_massZ],
        "AK4_OF" : [get_nMuon(1, coll="MuonGood"), get_nElectron(1, coll="ElectronGood"), get_nObj_min(4 , coll="CleanJet"),  get_nObj_eq(0, coll="CleanFatJet"), Vjet_massZ, dilepton_massZ],
        
        # DY
        "DoubleEle_AK8_AK4_sideL_Zjets" : [get_nElectron(1, coll="ElectronGood"),  Zjet_sideL],
        "DoubleEle_AK8_AK4_sideR_Zjets" : [get_nElectron(1, coll="ElectronGood"),  Zjet_sideR],
        "DoubleMuon_AK8_AK4_sideL_Zjets" : [get_nMuon(1, coll="MuonGood"),  Zjet_sideL],
        "DoubleMuon_AK8_AK4_sideR_Zjets" : [get_nMuon(1, coll="MuonGood"),  Zjet_sideR],
        "DoubleLepton_AK8_AK4_sideL_Zjets" : [get_nObj_eq(1, coll="LeptonGood"),  Zjet_sideL],
        "DoubleLepton_AK8_AK4_sideR_Zjets" : [get_nObj_eq(1, coll="LeptonGood"),  Zjet_sideR]

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
                           ColOut("dimuon_candidate", ["mass", "pt", "eta", "deltaEta", "deltaPhi", "deltaR", "charge"]),
                           ColOut("events", ["nCleanJets", "nCleanFatJets", "nBJetGood", "nJet", "nFatJet"], flatten=False),
                           ColOut("CleanSubJet_pair", ["zg"], flatten=False),
                           ],
            "bycategory" : {},
            },            
    },   
)