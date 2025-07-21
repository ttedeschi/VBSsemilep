from custom_cut_functions import *
from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.cut_definition import Cut
from pocket_coffea.lib.cut_functions import get_nObj_min, get_HLTsel, get_nPVgood, goldenJson, eventFlags, get_nElectron, get_nMuon, get_nObj_eq
from pocket_coffea.parameters.cuts import passthrough
from pocket_coffea.parameters.histograms import *
import workflowVBS
from workflowVBS import VBS_WV_Processor
from skim_dictionary import skim_dict
from pocket_coffea.lib.weights.common import common_weights
import cloudpickle
import custom_cut_functions 

cloudpickle.register_pickle_by_value(workflowVBS)
cloudpickle.register_pickle_by_value(custom_cut_functions)

import os
localdir = os.path.dirname(os.path.abspath(__file__))

# importing the dafult parameters for jec/jer...
from pocket_coffea.parameters import defaults
default_parameters = defaults.get_default_parameters()
defaults.register_configuration_dir("config_dir", localdir+"/parameters")

# categories...
# 1- electron (tight) + VBS jets
# 2- muon (tight) + VBS jets


parameters = defaults.merge_parameters_from_files(default_parameters,
                                    f"{localdir}/parameters/object_presel.yaml",
                                    update=True                                    
                                    )


# **************************************************************** #
# ************* MAIN CONFIGURATION ******************************* #
cfg = Configurator(
    parameters=parameters,
    datasets = {
        "jsons" : [f"{localdir}/datasets/VBS_mc_132X_Summer23wmLHEGS.json"],
        "filter" : {
            "samples" : ["ssWWLL"],
        }
        }, 
    workflow=VBS_WV_Processor,
     skim = [
            get_nPVgood(1), 
            get_nObj_min_or(
                [skim_dict["Wlep_V"]["Muon"]["ptMin"],skim_dict["Wlep_V"]["Electron"]["ptMin"]], 
                [skim_dict["Wlep_V"]["Muon"]["nMin"], skim_dict["Wlep_V"]["Electron"]["nMin"]],
                ["Muon", "Electron"]),
            ],
    preselections = [semileptonic_preselW, VBS_jets_presel],
    
    categories= {
        "baseline" : [passthrough],
        "SingleEle_plus_VBSjets" : [get_nElectron(1, coll="ElectronGood")],
        "SingleMuon_plus_VBSjets" : [get_nMuon(1, coll="MuonGood")],   
        "SingleEle_plus_VBSjets_oneAK8_twoAK4" : [get_nElectron(1, coll="ElectronGood"), get_nObj_eq(1, coll="CleanFatJets"), get_nObj_eq(2, coll="CleanJets")],
        "SingleEle_plus_VBSjets_fourAK4" : [get_nElectron(1, coll="ElectronGood"), get_nObj_eq(4 , coll="CleanJets")],
        "SingleMuon_plus_VBSjets_oneAK8_twoAK4" : [get_nMuon(1, coll="MuonGood"), get_nObj_eq(1, coll="CleanFatJets"), get_nObj_eq(2, coll="CleanJets")],
        "SingleMuon_plus_VBSjets_fourAK4" : [get_nMuon(1, coll="MuonGood"), get_nObj_eq(4 , coll="CleanJets")],
    },    
    
    weights_classes = common_weights,
    
    
    # weights/variations tutti da controllare in base a quelli disponibili per 2023
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


   variables = {
        **muon_hists(coll="MuonGood", pos=0),
        **count_hist(name="nElectronGood", coll="ElectronGood",bins=3, start=0, stop=3),
        **count_hist(name="nMuonGood", coll="MuonGood",bins=3, start=0, stop=3),
        **count_hist(name="nJets", coll="CleanJets",bins=8, start=0, stop=8),
        **count_hist(name="nBJets", coll="CleanFatJets",bins=8, start=0, stop=8),
    }

)