import ROOT
import math
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import *

class unpacking(Module):
    def __init__(self):
        pass
    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("Top_pt",  "F") 
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        muons = Collection(event, "Muon")
        jets = Collection(event, "Jet")
        goodMu = []
        goodJet = []
        top = ROOT.TLorentzVector()

        goodMu = list(filter(lambda x : x.pt>10, muons))
        goodJet = list(filter(lambda x : x.btagDeepFlavB > 0.5 , jets))

        if (len(goodMu)>0 and len(goodJet)>0):
            for k in goodMu:
                for j in goodJet:
                    top = j.p4() + k.p4()
                    self.out.fillBranch("Top_pt", top.Pt())

        return True

