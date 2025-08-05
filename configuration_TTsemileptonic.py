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
from pocket_coffea.lib.columns_manager import ColOut

cloudpickle.register_pickle_by_value(workflowVBS)
cloudpickle.register_pickle_by_value(custom_cut_functions)

import os
localdir = os.path.dirname(os.path.abspath(__file__))

# importing the dafult parameters for jec/jer...
from pocket_coffea.parameters import defaults
default_parameters = defaults.get_default_parameters()
defaults.register_configuration_dir("config_dir", localdir+"/parameters")

# ****************************************************************************** # 
# single-top with at least one lepton.







parameters = defaults.merge_parameters_from_files(default_parameters,
                                    f"{localdir}/parameters/object_presel.yaml",
                                    f"{localdir}/parameters/btagging.yaml",
                                    update=True                                    
                                    )


# **************************************************************** #
# ************* MAIN CONFIGURATION ******************************* #
cfg = Configurator(
    parameters=parameters,
    datasets = {
        "jsons" : ["/afs/cern.ch/user/l/ldellape/VBSsemilep/datasets/TTtoLNu2Q_HT-500_NJet-9_Hdamp-158_TuneCP5_13p6TeV_powheg-pythia8.json"],
        "filter" : {
            "samples" : ["TTtoLNu2Q_HT-500_NJet-9_Hdamp-158_TuneCP5_13p6TeV_powheg-pythia8"],
        }
        }, # here json files with the list of files
    workflow=VBS_WV_Processor,
     skim = [
            get_nPVgood(1), 
            get_nObj_min_or(
                [skim_dict["ttsemilep"]["Muon"]["ptMin"],skim_dict["ttsemilep"]["Electron"]["ptMin"]], 
                [skim_dict["ttsemilep"]["Muon"]["nMin"], skim_dict["ttsemilep"]["Electron"]["nMin"]],
                ["Muon", "Electron"]),
            ],
    preselections = [semileptonic_preselW, VBS_jets_presel, Bjets_presel],
    
    categories= {
        "baseline" : [passthrough],
        "SingleEle_plus_VBSjets_oneAK8_at_least_doubleB" : [get_nElectron(1, coll="ElectronGood"), get_nObj_min(2 , coll="BJetGood")],
        "SingleMuon_plus_VBSjets_oneAK8_at_least_doubleB" : [get_nMuon(1, coll="MuonGood"), get_nObj_min(2 , coll="BJetGood")],
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
    variables = {},
    workflow_options = {
        "dump_columns_as_arrays_per_chunk": "root://eosuser.cern.ch//eos/user/l/ldellape/VBS/parquet/"
    },
    columns = {
        "common" : {
            "inclusive" : [ColOut("ElectronGood", ["pt","eta", "phi"], flatten=False),
                           ColOut("MuonGood", ["pt", "eta", "phi"], flatten=False),
                           ColOut("CleanFatJet", ["pt", "eta", "phi", "tau1", "tau2", "tau21", "msoftdrop", "mass"], flatten=False),
                           ColOut("CleanJet", ["pt", "eta", "phi"], flatten=False),
                           ColOut("CleanJet_noVBS", ["pt", "eta", "phi"], flatten=False),
                           ColOut("VBS_dijet_system", ["mass", "pt", "deltaEta"], flatten=False),        
                           ColOut("BJetGood", ["mass", "pt", "eta", "phi"], flatten=False),
                           ],
            "bycategory" : {},
            },            
    },


)