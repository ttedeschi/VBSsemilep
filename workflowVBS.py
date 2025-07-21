import awkward as ak 
import numpy as np
from pocket_coffea.workflows.base import BaseProcessorABC
from pocket_coffea.utils.configurator import Configurator
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
        if self._isMC:
            self.events["Jet"] = jet_correction_correctionlib(self.events, "Jet", "AK4PFPuppi", "2023_Summer23", 'Summer23Prompt23_V2_MC') 
            self.events["FatJet"] = jet_correction_correctionlib(self.events, "FatJet", "AK8PFPuppi" , "2023_Summer23", 'Summer23Prompt23_V2_MC')
        
        # mask the electrons and muons
        self.events["MuonGood"] = lepton_selection(self.events, "Muon", self.params)
        self.events["ElectronGood"] = lepton_selection(self.events, "Electron", self.params)
        tight_muons = self.events["MuonGood"]
        
        
        print("**************")
        print(f" mvaIso : {self.events.Electron.mvaIso}")
        print(f" mvaIso_WP80 : {self.events.Electron.mvaIso_WP90}")
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
        sorted_indices = ak.argsort(self.events.CleanJet.pt, axis=1, ascending=False)
        print(f" ak4 before sorting: {self.events.CleanJet.pt}")
        sorted_CleanJet = self.events.CleanJet[sorted_indices]
        self.events["VBSjetsCandidate"] = sorted_CleanJet[:, :2]
        self.events["VBS_dijet_system"] = get_dijet(self.events["VBSjetsCandidate"], taggerVars=False)
        self.events["VjetsCandidate"] = sorted_CleanJet[:, 2:4]
        self.events["dijet_V_candidate"] = get_dijet(self.events["VjetsCandidate"], taggerVars=False)
        
        
        ''' non ci sono per 2023
        # now MET, checking for available corrections, if yes overwrite the collections
        MET_pt_corrected, MET_phi_corrected = met_xy_correction(self.params, self.events, self._year, self._era)
        self.events["MET"] = ak.with_field(self.events.MET, MET_pt_corrected, "pt")
        self.events["MET"] = ak.with_field(self.events.MET, MET_phi_corrected, "phi")
        '''
        
        
    def count_objects(self, variation):
        self.events["nMuonGood"] = ak.num(self.events.MuonGood)
        self.events["nElectronGood"] = ak.num(self.events.ElectronGood)
        self.events["nLeptonGood"] = ak.num(self.events.LeptonGood)
        self.events["nCleanFatJets"] = ak.num(self.events.CleanFatJet)
        self.events["nCleanJets"] = ak.num(self.events.CleanJet)
        print("**************************")
        print(f" n. good leptons: {self.events.nLeptonGood}")
        print(f" n. good muons: {self.events.nMuonGood}")
        print(f" n. good electrons: {self.events.nElectronGood}")
        print(f" n. clean jets {self.events.nCleanJets}")
        print(f" n. clean fatjets {self.events.nCleanFatJets}")