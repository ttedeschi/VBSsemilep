import awkward as ak 
import numpy as np
from pocket_coffea.workflows.base import BaseProcessorABC
from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.deltaR_matching import metric_eta, metric_phi
import vector

from pocket_coffea.lib.hist_manager import Axis
from correction_lib_jet import jet_correction_correctionlib
from pocket_coffea.lib.objects import (
    jet_correction,
    lepton_selection,
    jet_selection,
    btagging,
    get_dilepton,
    get_dijet,
    met_xy_correction
)

class VBS_WV_Processor(BaseProcessorABC):
    def __init__(self, cfg: Configurator):
        super().__init__(cfg)
        
    def apply_object_preselection(self, variation):
        # here the functions from pocket_coffea.lib.objects must be called
        # passing the yaml file defined in ./parameters directory
        if self._isMC and "2023" in self._year:
            self.events["Jet"] = jet_correction_correctionlib(self.events, "Jet", "AK4PFPuppi", "2023_Summer23", 'Summer23Prompt23_V2_MC') 
            self.events["FatJet"] = jet_correction_correctionlib(self.events, "FatJet", "AK8PFPuppi" , "2023_Summer23", 'Summer23Prompt23_V2_MC')
       
        
        # mask the electrons and muons
        self.events["MuonGood"] = lepton_selection(self.events, "Muon", self.params)
        self.events["ElectronGood"] = lepton_selection(self.events, "Electron", self.params)
        tight_muons = self.events["MuonGood"]
        
        
        tight_electron = self.events["ElectronGood"]
        self.events["LeptonGood"] = ak.concatenate((self.events.MuonGood, self.events.ElectronGood), axis=1)
        tight_leptons = self.events["LeptonGood"]
        
        
        # selection on jet and fatjet...checking leptons not overlapping with jets/fatjets
        self.events["CleanFatJet"], self.CleanFatJetMask = jet_selection(
            self.events, "FatJet", self.params, self._year,  leptons_collection="LeptonGood"
        )
        self.events["CleanJet"], self.CleanJetMask = jet_selection(
            self.events, "Jet", self.params, self._year,  leptons_collection="LeptonGood"
        )
        print("*****************************************************************************************")
        
        # create collection of the two ak4 jets that return the greatest mass and remove them from the 
        # ak4 list in a new CleanJet_noVBS colleciton
        self.VBS_pair_candidate()
        self.V_pair_candidate()
        self.Vlep_transverse()
        
        self.events["CleanFatJet"] = ak.with_field(
        self.events["CleanFatJet"],
        self.events.CleanFatJet.tau2 / self.events.CleanFatJet.tau1,
        "tau21"
        )      
        self.events["CleanFatJet"] = ak.with_field(
            self.events["CleanFatJet"],
            self.events.CleanFatJet.tau3/self.events.CleanFatJet.tau2,
            "tau32"
        )
        
        bjets_tagged = btagging(
            self.events["CleanJet"],
            self.params.btagging.working_point[self._year], 
            wp = self.params.object_preselection.Jet["btag"]["wp"],
        )
        self.events["BJetGood"] = bjets_tagged[abs(bjets_tagged.eta) < 2.5]

        
        ''' non ci sono per 2023
        # now MET, checking for available corrections, if yes overwrite the collections
        MET_pt_corrected, MET_phi_corrected = met_xy_correction(self.params, self.events, self._year, self._era)
        self.events["MET"] = ak.with_field(self.events.MET, MET_pt_corrected, "pt")
        self.events["MET"] = ak.with_field(self.events.MET, MET_phi_corrected, "phi")
        '''
        
        
    # make the pairs with the ak4 jets collection, take the pair with the greates m_jj as VBS dijet candidate
    def VBS_pair_candidate(self):
        ak4_cleanjets = self.events["CleanJet"]
        ak4_pairs = ak.combinations(ak4_cleanjets, 2, fields=["jet1", "jet2"])
        idx_pairs = ak.combinations(ak.local_index(ak4_cleanjets, axis=1), 2, fields=["idx1", "idx2"])
        dijet_mass = (ak4_pairs["jet1"] + ak4_pairs["jet2"]).mass
        sort_by_mass = ak.argsort(dijet_mass, axis=1, ascending=False)
        sorted_ak4_pair = ak4_pairs[sort_by_mass]
        sorted_ak4_idx = idx_pairs[sort_by_mass]
        best_pair = ak.firsts(sorted_ak4_pair)      
        best_idx_pair = ak.firsts(sorted_ak4_idx)    
        jet1_pt = best_pair["jet1"].pt
        jet2_pt = best_pair["jet2"].pt
        jet2_higher_pt = jet2_pt > jet1_pt
        jet1 = ak.where(jet2_higher_pt, best_pair["jet2"], best_pair["jet1"])
        jet2 = ak.where(jet2_higher_pt, best_pair["jet1"], best_pair["jet2"])
        idx1 = ak.where(jet2_higher_pt, best_idx_pair["idx2"], best_idx_pair["idx1"])
        idx2 = ak.where(jet2_higher_pt, best_idx_pair["idx1"], best_idx_pair["idx2"])
        dijet = best_pair["jet1"] + best_pair["jet2"]
        vbs_dijet = ak.zip({
            "mass": dijet.mass,
            "pt": dijet.pt,
            "deltaEta": abs(best_pair["jet1"].eta - best_pair["jet2"].eta),
            "deltaPhi": best_pair["jet1"].phi - best_pair["jet2"].phi,
            "pt1" : best_pair["jet1"].pt,
            "pt2" : best_pair["jet2"].pt,
            "idx1": best_idx_pair["idx1"],
            "idx2": best_idx_pair["idx2"],
        })
        self.events["VBS_dijet_system"] = vbs_dijet
        
    
    
    # create the pair of jets candidate from W/Z boson, removing the 2 jets already in the "VBS_dijet_system".
    # if more than two ak4 jets remains, take the pair with the sd mass closer to that of W
    def V_pair_candidate(self):
        ak4_cleanjets = self.events["CleanJet"]
        mask = ~(
            (ak.local_index(ak4_cleanjets, axis=1) == self.events["VBS_dijet_system"].idx1) |
            (ak.local_index(ak4_cleanjets, axis=1) == self.events["VBS_dijet_system"].idx2)
        )
        self.events["CleanJet_noVBS"] = self.events["CleanJet"][mask]
        clean_jets_no_vbs = self.events["CleanJet_noVBS"]

        # Count jets per event
        njets_remaining = ak.num(clean_jets_no_vbs)

        mask_eq2 = njets_remaining == 2
        mask_gt2 = njets_remaining > 2
        
        
        # only two ak4 remaining
        ak4_eq2 = ak.zip({
            "mass" : (clean_jets_no_vbs[mask_eq2][:,0] + clean_jets_no_vbs[mask_eq2][:,1]).mass,
            "pt" : (clean_jets_no_vbs[mask_eq2][:,0] + clean_jets_no_vbs[mask_eq2][:,1]).pt,
            "idx1" : ak.local_index(clean_jets_no_vbs[mask_eq2], axis=1)[:,0],
            "idx2" : ak.local_index(clean_jets_no_vbs[mask_eq2], axis=1)[:,1],
        }
        )
        
        # more than two ak4 remaining        
        jets_gt2 = clean_jets_no_vbs[mask_gt2]
        pair_gt2 = ak.combinations(jets_gt2, 2, fields=["jet1", "jet2"])
        idx_gt2 = ak.combinations(ak.local_index(jets_gt2, axis=1), 2, fields=["idx1", "idx2"])
        mass_gt2 = (pair_gt2["jet1"] + pair_gt2["jet2"]).mass
        diff80 = abs(mass_gt2 - 80)
        sort_idx = ak.argsort(diff80, axis=1, ascending=True)
        sorted_pair = pair_gt2[sort_idx]
        sorted_idx = idx_gt2[sort_idx]
        best_pair_gt2 = ak.firsts(sorted_pair)
        best_idx_gt2 = ak.firsts(sorted_idx)
        dijet_gt2 = best_pair_gt2["jet1"] + best_pair_gt2["jet2"]
        result_gt2 = ak.zip({
            "mass": dijet_gt2.mass,
            "pt": dijet_gt2.pt,
            "deltaEta": abs(best_pair_gt2["jet1"].eta - best_pair_gt2["jet2"].eta),
            "deltaPhi": best_pair_gt2["jet1"].phi - best_pair_gt2["jet2"].phi,
            "idx1": best_idx_gt2["idx1"],
            "idx2": best_idx_gt2["idx2"],
        })
        n_events = len(self.events)
        full_result = [None] * n_events 

        ak4_eq2_idx = ak.where(mask_eq2)[0]  
        for i, val in zip(ak4_eq2_idx, ak4_eq2):
            if val.mass is not None:
                full_result[i] = val

        result_gt2_idx = ak.where(mask_gt2)[0]  
        for i, val in zip(result_gt2_idx, result_gt2):
            if val.mass is not None:
                full_result[i] = val

        self.events["V_dijet_candidate"] = ak.Array(full_result)
        
        
    # prepare the transverse mass of electron/muon + MET
    def Vlep_transverse(self):
        self.events["MT_lep_miss"] = np.sqrt(
            2 * self.events.LeptonGood.pt * self.events.MET.pt *
            (1 - np.cos(self.events.LeptonGood.phi - self.events.MET.phi))
        )
        self.events["MT_ele_miss"] = np.sqrt(
            2 * self.events.ElectronGood.pt * self.events.MET.pt *
            (1 - np.cos(self.events.ElectronGood.phi - self.events.MET.phi))
        )
        self.events["MT_mu_miss"] = np.sqrt(
            2 * self.events.MuonGood.pt * self.events.MET.pt *
            (1 - np.cos(self.events.MuonGood.phi - self.events.MET.phi))
        )
        
        
        
    def count_objects(self, variation):
        self.events["nMuonGood"] = ak.num(self.events.MuonGood)
        self.events["nElectronGood"] = ak.num(self.events.ElectronGood)
        self.events["nLeptonGood"] = ak.num(self.events.LeptonGood)
        self.events["nCleanFatJets"] = ak.num(self.events.CleanFatJet)
        self.events["nCleanJets"] = ak.num(self.events.CleanJet)
        self.events["nBJetGood"] = ak.num(self.events.BJetGood)
