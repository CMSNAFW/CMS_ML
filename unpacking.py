import ROOT
import math
import numpy as np
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import *
from PhysicsTools.NanoAODTools.postprocessing.skimtree_utils import *

class unpacking(Module):
    def __init__(self):
        pass
    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        """Branch variabili top senza nu"""
        self.out.branch("Top_pt","F", lenVar="top.Size()") 
        self.out.branch("Top_eta","F", lenVar="top.Size()")
        self.out.branch("Top_phi","F", lenVar="top.Size()")
        self.out.branch("Top_e","F", lenVar="top.Size()")
        self.out.branch("Top_M","F", lenVar="top.Size()")

        """Branch variabili top con nu"""

        self.out.branch("Top_nu_pt","F", lenVar="top.Size()") 
        self.out.branch("Top_nu_eta","F", lenVar="top.Size()")
        self.out.branch("Top_nu_phi","F", lenVar="top.Size()")
        self.out.branch("Top_nu_e","F", lenVar="top.Size()")
        self.out.branch("Top_nu_M","F", lenVar="top.Size()")

        """Branch indici goodJet e goodMu"""
        self.out.branch("Top_bjet_index","i", lenVar="top.Size()")
        self.out.branch("Top_mu_index","i", lenVar="top.Size()")

        """Muon and bjet unboosted (top frame)"""

        self.out.branch("Jet_unboosted_pt","F", lenVar="top.Size()") 
        self.out.branch("Jet_unboosted_eta","F", lenVar="top.Size()")
        self.out.branch("Jet_unboosted_phi","F", lenVar="top.Size()")
        self.out.branch("Jet_unboosted_e","F", lenVar="top.Size()")
        self.out.branch("Jet_unboosted_M","F", lenVar="top.Size()")

        self.out.branch("Muon_unboosted_pt","F", lenVar="top.Size()") 
        self.out.branch("Muon_unboosted_eta","F", lenVar="top.Size()")
        self.out.branch("Muon_unboosted_phi","F", lenVar="top.Size()")
        self.out.branch("Muon_unboosted_e","F", lenVar="top.Size()")
        self.out.branch("Muon_unboosted_M","F", lenVar="top.Size()")

        self.out.branch("Top_pt_rel","F", lenVar="top.Size()")
        self.out.branch("Is_dR_merg","B", lenVar="top.Size()")
        self.out.branch("Costheta","F", lenVar="top.Size()")

        self.out.branch("Top_High_Truth","F", lenVar="top.Size()")
        self.out.branch("Tau_High_Truth","F", lenVar="top.Size()")


    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        muons = Collection(event, "Muon")
        jets = Collection(event, "Jet")
        MET = Object(event, "MET")
        genpart = Collection(event, "GenPart")
        goodMu = []
        goodJet = []
        """Variabili top senza nu"""
        top_momentum=ROOT.TLorentzVector()
        top_pt =[]
        top_phi = []
        top_eta = []
        top_e = []
        top_M = []

        """Variabili top con nu"""
        top_nu_momentum_utils= TopUtilities()
        top_nu_pt =[]
        top_nu_phi = []
        top_nu_eta = []
        top_nu_e = []
        top_nu_M = []

        top_nu_momentum = ROOT.TLorentzVector()
        top_nu_momentum_neg = ROOT.TLorentzVector()
        IsmcNeg = False
        mcdR_lepjet = None
        costheta = []

        """variabili bjet e muon"""
        top_bjet_index = []
        top_mu_index = []

        muon_unboosted_momentum = ROOT.TLorentzVector()
        muon_unboosted_pt =[]
        muon_unboosted_phi = []
        muon_unboosted_eta = []
        muon_unboosted_e = []
        muon_unboosted_M = []

        jet_unboosted_momentum = ROOT.TLorentzVector()
        jet_unboosted_pt =[]
        jet_unboosted_phi = []
        jet_unboosted_eta = []
        jet_unboosted_e = []
        jet_unboosted_M = []

        top_pt_rel = []
        is_dR_merg = []

        top_high_truth = []
        tau_high_truth = []

        """"""

        bjets, nobjets = bjet_filter(jets, 'DeepCSV', 'L')
        
        goodMu = list(filter(lambda x : x.pt>10, muons))
        goodJet = list(filter(lambda x : x.pt>10 , bjets))

        if not (len(goodMu)>0 and len(goodJet)>0):
            
            return False
        
        if (len(goodMu)>0 and len(goodJet)>0):
            for k in goodMu:
                for j in goodJet:

                    is_muon_prompt = False
                    is_jet_true = False
                    is_muon_from_tau = False
                    
                    top_momentum = (j.p4() + k.p4())
                    top_pt.append(top_momentum.Pt())
                    top_phi.append(top_momentum.Phi())
                    top_eta.append(top_momentum.Eta())
                    top_e.append(top_momentum.E())
                    top_M.append(top_momentum.M())
 
                    top_nu_momentum, IsmcNeg, mcdR_lepjet = top_nu_momentum_utils.top4Momentum(k.p4(),j.p4(),MET.pt*math.cos(MET.phi),MET.pt*math.sin(MET.phi))

                    top_nu_pt.append(top_nu_momentum.Pt())
                    top_nu_phi.append(top_nu_momentum.Phi())
                    top_nu_eta.append(top_nu_momentum.Eta())
                    top_nu_e.append(top_nu_momentum.E())
                    top_nu_M.append(top_nu_momentum.M())  

                    top_bjet_index.append(j)
                    top_mu_index.append(k)  
                    top_pt_rel.append(((k.p4().Vect()).Cross(j.p4().Vect())).Mag()/((j.p4().Vect()).Mag())) 

                    """unboosting"""
                    top_nu_momentum_neg.SetPxPyPzE(-top_nu_momentum.Px(),-top_nu_momentum.Py(),-top_nu_momentum.Pz(),top_nu_momentum.E())

                    muon_unboosted_momentum = k.p4()
                    muon_unboosted_momentum.Boost(top_nu_momentum_neg.BoostVector())
                    muon_unboosted_pt.append(muon_unboosted_momentum.Pt())
                    muon_unboosted_phi.append(muon_unboosted_momentum.Phi())
                    muon_unboosted_eta.append(muon_unboosted_momentum.Eta())
                    muon_unboosted_e.append(muon_unboosted_momentum.E())
                    muon_unboosted_M.append(muon_unboosted_momentum.M())

                    jet_unboosted_momentum = j.p4()
                    jet_unboosted_momentum.Boost(top_nu_momentum_neg.BoostVector())

                    jet_unboosted_pt.append(jet_unboosted_momentum.Pt())
                    jet_unboosted_phi.append(jet_unboosted_momentum.Phi())
                    jet_unboosted_eta.append(jet_unboosted_momentum.Eta())
                    jet_unboosted_e.append(jet_unboosted_momentum.E())
                    jet_unboosted_M.append(jet_unboosted_momentum.M())

                    if deltaR(j.p4().Eta(),j.p4().Phi(),k.p4().Eta(),k.p4().Phi()) <0.4 :
                        is_dR_merg.append(False)
                    else:
                        is_dR_merg.append(True)

                    costheta.append(top_nu_momentum_utils.costhetapol(k.p4(),j.p4(),top_nu_momentum))

                    if (k.genPartFlav== 1) : 
                        for gen in genpart:
                            if (gen == k.genPartIdx):
                                if(genpart.pdgId[gen.genPartIdxMother]==24):
                                    is_muon_prompt = True

                    if (k.genPartFlav== 15) : 
                        is_muon_from_tau = True

                    if ((j.partonFlavour)==5):
                        is_jet_true = True

                    if ((is_jet_true * is_muon_prompt)== True and (j.partonFlavour*k.charge)<0.) :
                        top_high_truth.append(5)

                    if (is_jet_true == False and is_muon_prompt== False):
                        top_high_truth.append(0)

                    if (is_jet_true == True and is_muon_prompt== False):
                        top_high_truth.append(1)


                    
                    if (is_jet_true == False and is_muon_prompt== True):
                        top_high_truth.append(3)


                    if ((is_jet_true * is_muon_prompt)== True and (j.partonFlavour*k.charge)>0.) :
                        top_high_truth.append(4)

                    tau_high_truth.append(is_muon_from_tau)

                    
        self.out.fillBranch("Top_pt", top_pt)
        self.out.fillBranch("Top_phi", top_phi)
        self.out.fillBranch("Top_eta", top_eta)
        self.out.fillBranch("Top_e", top_e)
        self.out.fillBranch("Top_M", top_M)

        self.out.fillBranch("Top_nu_pt", top_nu_pt)
        self.out.fillBranch("Top_nu_phi", top_nu_phi)
        self.out.fillBranch("Top_nu_eta", top_nu_eta)
        self.out.fillBranch("Top_nu_e", top_nu_e)
        self.out.fillBranch("Top_nu_M", top_nu_M)

        #self.out.fillBranch("Top_bjet_index", top_bjet_index)
        #self.out.fillBranch("Top_mu_index", top_mu_index)
        self.out.fillBranch("Top_pt_rel", top_pt_rel)

        self.out.fillBranch("Jet_unboosted_pt",jet_unboosted_pt) 
        self.out.fillBranch("Jet_unboosted_eta",jet_unboosted_eta)
        self.out.fillBranch("Jet_unboosted_phi",jet_unboosted_phi)
        self.out.fillBranch("Jet_unboosted_e",jet_unboosted_e)
        self.out.fillBranch("Jet_unboosted_M",jet_unboosted_M)

        self.out.fillBranch("Muon_unboosted_pt",muon_unboosted_pt) 
        self.out.fillBranch("Muon_unboosted_eta",muon_unboosted_eta)
        self.out.fillBranch("Muon_unboosted_phi",muon_unboosted_phi)
        self.out.fillBranch("Muon_unboosted_e",muon_unboosted_e)
        self.out.fillBranch("Muon_unboosted_M",muon_unboosted_M)

        self.out.fillBranch("Is_dR_merg",is_dR_merg)
        self.out.fillBranch("Costheta", costheta)

        self.out.fillBranch("Top_High_Truth",top_high_truth)
        self.out.fillBranch("Tau_High_Truth",tau_high_truth)


        return True

