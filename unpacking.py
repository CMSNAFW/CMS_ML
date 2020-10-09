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


    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        muons = Collection(event, "Muon")
        jets = Collection(event, "Jet")
        MET = Object(event, "MET")
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
        IsmcNeg = False
        mcdR_lepjet = None


        bjets, nobjets = bjet_filter(jets, 'DeepFlv', 'L')
        
        goodMu = list(filter(lambda x : x.pt>10, muons))
        goodJet = list(filter(lambda x : x.pt>10 , bjets))

        if not (len(goodMu)>0 and len(goodJet)>0):
            
            return False
        
        if (len(goodMu)>0 and len(goodJet)>0):
            for k in goodMu:
                for j in goodJet:
                    
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

        return True

