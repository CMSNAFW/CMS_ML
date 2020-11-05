import ROOT
import math
import numpy as np
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import *
from PhysicsTools.NanoAODTools.postprocessing.skimtree_utils import *

class unpacking_vers2(Module):
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
        self.out.branch("Top_bjet_index","F", lenVar="top.Size()")
        self.out.branch("Top_mu_index","F", lenVar="top.Size()")
        self.out.branch("Top_el_index","F", lenVar = "top.Size()")

        """Muon and bjet unboosted (top frame)"""

        self.out.branch("Jet_unboosted_pt","F", lenVar="top.Size()") 
        self.out.branch("Jet_unboosted_eta","F", lenVar="top.Size()")
        self.out.branch("Jet_unboosted_phi","F", lenVar="top.Size()")
        self.out.branch("Jet_unboosted_e","F", lenVar="top.Size()")
        self.out.branch("Jet_unboosted_M","F", lenVar="top.Size()")
        self.out.branch("Jet_has_promptLep","F", lenVar="top.Size()")

        self.out.branch("Lep_unboosted_pt","F", lenVar="top.Size()") 
        self.out.branch("Lep_unboosted_eta","F", lenVar="top.Size()")
        self.out.branch("Lep_unboosted_phi","F", lenVar="top.Size()")
        self.out.branch("Lep_unboosted_e","F", lenVar="top.Size()")
        self.out.branch("Lep_unboosted_M","F", lenVar="top.Size()")

        self.out.branch("Top_pt_rel","F", lenVar="top.Size()")
        self.out.branch("Is_dR_merg","F", lenVar="top.Size()")
        self.out.branch("Costheta","F", lenVar="top.Size()")
        self.out.branch("Top_dR", "F",lenVar ="top.Size()")

        self.out.branch("Top_High_Truth","F", lenVar="top.Size()")
        self.out.branch("Tau_High_Truth","B", lenVar="top.Size()")


    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        muons = Collection(event, "Muon")
        jets = Collection(event, "Jet")
        electrons = Collection(event, "Electron")
        MET = Object(event, "MET")
        genpart = Collection(event, "GenPart")
        LHE = Collection(event, "LHEPart")
    
        
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
        top_el_index = []

        lep_unboosted_momentum = ROOT.TLorentzVector()
        lep_unboosted_pt =[]
        lep_unboosted_phi = []
        lep_unboosted_eta = []
        lep_unboosted_e = []
        lep_unboosted_M = []

        jet_unboosted_momentum = ROOT.TLorentzVector()
        jet_unboosted_pt =[]
        jet_unboosted_phi = []
        jet_unboosted_eta = []
        jet_unboosted_e = []
        jet_unboosted_M = []

        top_pt_rel = []
        is_dR_merg = []
        top_dR = []

        top_high_truth = []
        tau_high_truth = []
        jet_has_pL = []
       

        """"""

        bjets, nobjets = bjet_filter(jets, 'DeepCSV', 'L')
        
        goodMu = list(filter(lambda x : x.pt>10, muons))
        goodJet = list(filter(lambda x :  x.pt>10, jets)) #noDeePCSV
        goodEl = list(filter(lambda x : x.pt>10, electrons))

        allMu = list(filter(lambda x : x.pt>0, muons))
        allJet = list(filter(lambda x :  x.pt>0 , jets))
        allEl = list(filter(lambda x : x.pt>0, electrons))

        if not ((len(goodMu)>0 or len(goodEl)) and len(goodJet)>0 ):
            
            return False
        
        else:


            for j in allJet:  

               


                if j in goodJet:

                    is_jet_true = False
                    jet_has_promptLep = False 
                    
                    for gen in LHE:
                            if (deltaR(j.p4().Eta(),j.p4().Phi(),gen.p4().Eta(),gen.p4().Phi()) <0.4 and abs(gen.pdgId)==5 ):
                                #if(abs(genpart[gen.genPartIdxMother].pdgId)==6 ):
                                    is_jet_true = True

                       
                    for gen in LHE:
                        if((abs(gen.pdgId) == 11 or abs(gen.pdgId)==13 or abs(gen.pdgId)== 15)):
                            #if(abs(genpart[gen.genPartIdxMother].pdgId)==24 and genpart[gen.genPartIdxMother].genPartIdxMother>-1): 
                                #if (abs(genpart[genpart[gen.genPartIdxMother].genPartIdxMother].pdgId)==6 ):
                                    if deltaR(j.p4().Eta(),j.p4().Phi(),gen.p4().Eta(),gen.p4().Phi()) <0.4 :
                                        jet_has_promptLep = True 
                    jet_has_pL.append(jet_has_promptLep)
                    
                    for m in allMu:

                        if m in goodMu:
                            is_muon_prompt = False
                    
                            is_muon_from_tau = False
                           
                            top_bjet_index.append(allJet.index(j))
                            top_mu_index.append(allMu.index(m))
                            top_el_index.append(-1)


                            if deltaR(j.p4().Eta(),j.p4().Phi(),m.p4().Eta(),m.p4().Phi()) <=0.4 :
                                is_dR_merg.append(1)


                            elif( deltaR(j.p4().Eta(),j.p4().Phi(),m.p4().Eta(),m.p4().Phi()) <=2 and deltaR(j.p4().Eta(),j.p4().Phi(),m.p4().Eta(),m.p4().Phi()) >0.4):
                                is_dR_merg.append(0)

                            else: 
                                is_dR_merg.append(-1)
                                


          
                            top_momentum = (j.p4() + m.p4())
                            top_pt.append(top_momentum.Pt())
                            top_phi.append(top_momentum.Phi())
                            top_eta.append(top_momentum.Eta())
                            top_e.append(top_momentum.E())
                            top_M.append(top_momentum.M())

                            

         
                            top_nu_momentum, IsmcNeg, mcdR_lepjet = top_nu_momentum_utils.top4Momentum(m.p4(),j.p4(),MET.pt*math.cos(MET.phi),MET.pt*math.sin(MET.phi))

                            top_nu_pt.append(top_nu_momentum.Pt())
                            top_nu_phi.append(top_nu_momentum.Phi())
                            top_nu_eta.append(top_nu_momentum.Eta())
                            top_nu_e.append(top_nu_momentum.E())
                            top_nu_M.append(top_nu_momentum.M())  

          
                            top_pt_rel.append(((m.p4().Vect()).Cross(j.p4().Vect())).Mag()/((j.p4().Vect()).Mag())) 

                            """unboosting"""
                            top_nu_momentum_neg.SetPxPyPzE(-top_nu_momentum.Px(),-top_nu_momentum.Py(),-top_nu_momentum.Pz(),top_nu_momentum.E())

                            lep_unboosted_momentum = m.p4()
                            lep_unboosted_momentum.Boost(top_nu_momentum_neg.BoostVector())
                            lep_unboosted_pt.append(lep_unboosted_momentum.Pt())
                            lep_unboosted_phi.append(lep_unboosted_momentum.Phi())
                            lep_unboosted_eta.append(lep_unboosted_momentum.Eta())
                            lep_unboosted_e.append(lep_unboosted_momentum.E())
                            lep_unboosted_M.append(lep_unboosted_momentum.M())

                            jet_unboosted_momentum = j.p4()
                            jet_unboosted_momentum.Boost(top_nu_momentum_neg.BoostVector())

                            jet_unboosted_pt.append(jet_unboosted_momentum.Pt())
                            jet_unboosted_phi.append(jet_unboosted_momentum.Phi())
                            jet_unboosted_eta.append(jet_unboosted_momentum.Eta())
                            jet_unboosted_e.append(jet_unboosted_momentum.E())
                            jet_unboosted_M.append(jet_unboosted_momentum.M())



                            top_dR.append(deltaR(j.p4().Eta(),j.p4().Phi(),m.p4().Eta(),m.p4().Phi()))

                            costheta.append(top_nu_momentum_utils.costhetapol(m.p4(),j.p4(),top_nu_momentum))

                            for gen in LHE:
                                if ( abs(gen.pdgId)==13 ):
                                    if(deltaR(m.p4().Eta(),m.p4().Phi(),gen.p4().Eta(),gen.p4().Phi()) <0.4):
                                #if(abs(genpart[gen.genPartIdxMother].pdgId)==6 ):
                                     is_muon_prompt = True
                            
                            for gen in LHE:
                                if (deltaR(m.p4().Eta(),m.p4().Phi(),gen.p4().Eta(),gen.p4().Phi()) <0.4 and abs(gen.pdgId)==15 ):
                                #if(abs(genpart[gen.genPartIdxMother].pdgId)==6 ):
                                    is_muon_from_tau = True
                                        

                 
                                        

                            tau_high_truth.append(is_muon_from_tau)

                            if ((is_jet_true * is_muon_prompt)== True and (j.partonFlavour*m.charge)>0.) :
                                top_high_truth.append(0)

                            if (is_jet_true == False and is_muon_prompt== False):
                                top_high_truth.append(1)

                            if (is_jet_true == True and is_muon_prompt== False and jet_has_promptLep == True):
                                top_high_truth.append(2)

                            if (is_jet_true == True and is_muon_prompt== False and jet_has_promptLep == False):
                                top_high_truth.append(3)
                            
                            if (is_jet_true == False and is_muon_prompt== True):
                                top_high_truth.append(4)


                            if ((is_jet_true * is_muon_prompt)== True and (j.partonFlavour*m.charge)<0.) :
                                top_high_truth.append(5)


                    for e in allEl:

                        if e in goodEl:
                            is_el_prompt = False
                    
                            is_el_from_tau = False
                           
                            top_bjet_index.append(allJet.index(j))
                            top_mu_index.append(-1)
                            top_el_index.append(allEl.index(e))


                            if deltaR(j.p4().Eta(),j.p4().Phi(),e.p4().Eta(),e.p4().Phi()) <=0.4 :
                                is_dR_merg.append(1)


                            elif( deltaR(j.p4().Eta(),j.p4().Phi(),e.p4().Eta(),e.p4().Phi()) <=2 and deltaR(j.p4().Eta(),j.p4().Phi(),e.p4().Eta(),e.p4().Phi()) >0.4):
                                is_dR_merg.append(0)

                            else: 
                                is_dR_merg.append(-1)
                                


          
                            top_momentum = (j.p4() + e.p4())
                            top_pt.append(top_momentum.Pt())
                            top_phi.append(top_momentum.Phi())
                            top_eta.append(top_momentum.Eta())
                            top_e.append(top_momentum.E())
                            top_M.append(top_momentum.M())

                           

         
                            top_nu_momentum, IsmcNeg, mcdR_lepjet = top_nu_momentum_utils.top4Momentum(e.p4(),j.p4(),MET.pt*math.cos(MET.phi),MET.pt*math.sin(MET.phi))

                            top_nu_pt.append(top_nu_momentum.Pt())
                            top_nu_phi.append(top_nu_momentum.Phi())
                            top_nu_eta.append(top_nu_momentum.Eta())
                            top_nu_e.append(top_nu_momentum.E())
                            top_nu_M.append(top_nu_momentum.M())  

          
                            top_pt_rel.append(((e.p4().Vect()).Cross(j.p4().Vect())).Mag()/((j.p4().Vect()).Mag())) 

                            """unboosting"""
                            top_nu_momentum_neg.SetPxPyPzE(-top_nu_momentum.Px(),-top_nu_momentum.Py(),-top_nu_momentum.Pz(),top_nu_momentum.E())

                            lep_unboosted_momentum = e.p4()
                            lep_unboosted_momentum.Boost(top_nu_momentum_neg.BoostVector())
                            lep_unboosted_pt.append(lep_unboosted_momentum.Pt())
                            lep_unboosted_phi.append(lep_unboosted_momentum.Phi())
                            lep_unboosted_eta.append(lep_unboosted_momentum.Eta())
                            lep_unboosted_e.append(lep_unboosted_momentum.E())
                            lep_unboosted_M.append(lep_unboosted_momentum.M())

                            jet_unboosted_momentum = j.p4()
                            jet_unboosted_momentum.Boost(top_nu_momentum_neg.BoostVector())

                            jet_unboosted_pt.append(jet_unboosted_momentum.Pt())
                            jet_unboosted_phi.append(jet_unboosted_momentum.Phi())
                            jet_unboosted_eta.append(jet_unboosted_momentum.Eta())
                            jet_unboosted_e.append(jet_unboosted_momentum.E())
                            jet_unboosted_M.append(jet_unboosted_momentum.M())



                            top_dR.append(deltaR(j.p4().Eta(),j.p4().Phi(),e.p4().Eta(),e.p4().Phi()))

                            costheta.append(top_nu_momentum_utils.costhetapol(e.p4(),j.p4(),top_nu_momentum))

 
                            for gen in LHE:
                                if (deltaR(e.p4().Eta(),e.p4().Phi(),gen.p4().Eta(),gen.p4().Phi()) <0.4 and abs(gen.pdgId)==11 ):
                                #if(abs(genpart[gen.genPartIdxMother].pdgId)==6 ):
                                    is_el_prompt = True
                            
                            for gen in LHE:
                                if (deltaR(e.p4().Eta(),e.p4().Phi(),gen.p4().Eta(),gen.p4().Phi()) <0.4 and abs(gen.pdgId)==15 ):
                                #if(abs(genpart[gen.genPartIdxMother].pdgId)==6 ):
                                    is_el_from_tau = True

                            tau_high_truth.append(is_el_from_tau)

                            if (is_jet_true== True  and is_el_prompt== True and (j.partonFlavour*e.charge)>0.) :
                                top_high_truth.append(0)

                            if (is_jet_true == False and is_el_prompt== False):
                                top_high_truth.append(1)

                            if (is_jet_true == True and is_el_prompt== False and jet_has_promptLep == True):
                                top_high_truth.append(2)

                            if (is_jet_true == True and is_el_prompt== False and jet_has_promptLep == False):
                                top_high_truth.append(3)
                            
                            if (is_jet_true == False and is_el_prompt== True):
                                top_high_truth.append(4)


                            if ((is_jet_true * is_el_prompt)== True and (j.partonFlavour*e.charge)<0.) :
                                top_high_truth.append(5)

                    
                    





                    
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

        self.out.fillBranch("Top_bjet_index", top_bjet_index)
        self.out.fillBranch("Top_mu_index", top_mu_index)
        self.out.fillBranch("Top_el_index", top_el_index)
        self.out.fillBranch("Top_pt_rel", top_pt_rel)

        self.out.fillBranch("Jet_unboosted_pt",jet_unboosted_pt) 
        self.out.fillBranch("Jet_unboosted_eta",jet_unboosted_eta)
        self.out.fillBranch("Jet_unboosted_phi",jet_unboosted_phi)
        self.out.fillBranch("Jet_unboosted_e",jet_unboosted_e)
        self.out.fillBranch("Jet_unboosted_M",jet_unboosted_M)
        self.out.fillBranch("Jet_has_promptLep",jet_has_pL)

        self.out.fillBranch("Lep_unboosted_pt",lep_unboosted_pt) 
        self.out.fillBranch("Lep_unboosted_eta",lep_unboosted_eta)
        self.out.fillBranch("Lep_unboosted_phi",lep_unboosted_phi)
        self.out.fillBranch("Lep_unboosted_e",lep_unboosted_e)
        self.out.fillBranch("Lep_unboosted_M",lep_unboosted_M)

        self.out.fillBranch("Is_dR_merg",is_dR_merg)
        self.out.fillBranch("Costheta", costheta)
        self.out.fillBranch("Top_dR", top_dR)

        self.out.fillBranch("Top_High_Truth",top_high_truth)
        self.out.fillBranch("Tau_High_Truth",tau_high_truth)


        return True

