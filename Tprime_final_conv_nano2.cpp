#include "TMath.h"
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <math.h>
#include <time.h>
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include <vector>
#include <assert.h>
#include <TMVA/Reader.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cassert>
#include <sstream>
#include <string>
#include "TFileCollection.h"
#include "THashList.h"
#include "TBenchmark.h"
#include "../src/DMTopVariables.h"
//#include "/afs/cern.ch/user/f/fcarneva/CMSSW_10_2_0/src/Tprime/TprimeAnalysis/src/DMTopVariables.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "interface/Weights.h"
#include "interface/MT2Utility.h"
#include "interface/mt2w_bisect.h"
#include "interface/mt2bl_bisect.h"
#include "interface/Mt2Com_bisect.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "interface/topTagging.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"


void Tprime_top1(string filename, string output, int f){
	
	 gStyle->SetOptStat("e");
  gStyle->SetStatX(0.899514);
  gStyle->SetStatY(0.901);
  gStyle->SetStatW(0.3); 
  gStyle->SetStatH(0.25);




  TFile *f1 =  TFile::Open((filename).c_str());
  TTree *tree = new TTree;
    tree = (TTree*) f1->Get("Events");
  int nentries = tree->GetEntries(); 
 

int size_max=20;
int size_max2=300;

 Float_t mu_pt[size_max]; tree->SetBranchAddress("Muon_pt",&mu_pt);
 Float_t mu_mass[size_max]; tree->SetBranchAddress("Muon_mass",&mu_mass);
 Float_t mu_phi[size_max]; tree->SetBranchAddress("Muon_phi",&mu_phi);
 Float_t mu_eta[size_max]; tree->SetBranchAddress("Muon_eta",&mu_eta);
 UInt_t musize = tree->SetBranchAddress("nMuon",&musize);
 Int_t mu_ch[size_max]; tree->SetBranchAddress("Muon_charge",&mu_ch);
 //Float_t mu_DB[size_max]; tree->SetBranchAddress("Muon_DB",&mu_DB);
  Float_t mu_Dxyerr[size_max]; tree->SetBranchAddress("Muon_dxyErr",&mu_Dxyerr);
   Float_t mu_Dz[size_max]; tree->SetBranchAddress("Muon_dz",&mu_Dz);
 Float_t mu_Dxy[size_max]; tree->SetBranchAddress("Muon_dxy",&mu_Dxy);
 Float_t mu_MiniIso[size_max]; tree->SetBranchAddress("Muon_miniPFRelIso_all",&mu_MiniIso);
 Float_t mu_Iso[size_max]; tree->SetBranchAddress("Muon_pfRelIso04_all",&mu_Iso);
  Float_t mu_isHigh[size_max]; tree->SetBranchAddress("Muon_highPtId",&mu_isHigh);
  Float_t mu_isTight[size_max]; tree->SetBranchAddress("Muon_tightId",&mu_isTight);
 Float_t mu_IsGlobal[size_max]; tree->SetBranchAddress("Muon_isGlobal",&mu_IsGlobal);
 Float_t mu_IsTracker[size_max]; tree->SetBranchAddress("Muon_isTracker",&mu_IsTracker);
 Float_t mu_NumberMatchedStations[size_max]; tree->SetBranchAddress("Muon_nStations",&mu_NumberMatchedStations);
UChar_t mu_match[size_max]; tree->SetBranchAddress("Muon_genPartFlav",&mu_match);
Int_t mu_gen_idx[size_max]; tree->SetBranchAddress("Muon_genPartIdx",&mu_gen_idx); 


float met_pt= tree->SetBranchAddress("MET_pt",&met_pt);
float met_phi= tree->SetBranchAddress("MET_phi",&met_phi);


 Float_t jet_pt[size_max]; tree->SetBranchAddress("Jet_pt",&jet_pt);
 Float_t jet_mass[size_max]; tree->SetBranchAddress("Jet_mass",&jet_mass);
 Float_t jet_phi[size_max]; tree->SetBranchAddress("Jet_phi",&jet_phi);
 Float_t jet_eta[size_max]; tree->SetBranchAddress("Jet_eta",&jet_eta);
 UInt_t jetsize = tree->SetBranchAddress("nJet",&jetsize);
 Int_t jet_fl[size_max]; tree->SetBranchAddress("Jet_hadronFlavour",&jet_fl);
Float_t jet_isDeep[size_max]; tree->SetBranchAddress("Jet_btagDeepFlavB",&jet_isDeep); 
Int_t Jet_partonFlavour[size_max]; tree->SetBranchAddress("Jet_partonFlavour",&Jet_partonFlavour);
Int_t Jet_genJetIdx[size_max]; tree->SetBranchAddress("Jet_genJetIdx",&Jet_genJetIdx);

 UInt_t gP_size= tree->SetBranchAddress("nGenPart",&gP_size);
 Float_t gP_Pt[size_max2]; tree->SetBranchAddress("GenPart_pt",&gP_Pt);
 Float_t gP_Phi[size_max2]; tree->SetBranchAddress("GenPart_phi",&gP_Phi);
 Float_t gP_Eta[size_max2]; tree->SetBranchAddress("GenPart_eta",&gP_Eta);
 Float_t gP_mass[size_max2]; tree->SetBranchAddress("GenPart_mass",&gP_mass);
 Int_t gP_flavor[size_max2]; tree->SetBranchAddress("GenPart_pdgId",&gP_flavor);
 Int_t gP_mom1I[size_max2]; tree->SetBranchAddress("GenPart_genPartIdxMother",&gP_mom1I);

 
Float_t FatJet_deepTag_H[size_max]; tree->SetBranchAddress("FatJet_deepTag_H",&FatJet_deepTag_H);
Float_t FatJet_deepTagMD_H4qvsQCD[size_max]; tree->SetBranchAddress("FatJet_deepTagMD_H4qvsQCD",&FatJet_deepTagMD_H4qvsQCD);
Float_t FatJet_deepTagMD_HbbvsQCD[size_max]; tree->SetBranchAddress("FatJet_deepTagMD_HbbvsQCD",&FatJet_deepTagMD_HbbvsQCD);
Float_t FatJet_deepTagMD_ZHbbvsQCD[size_max]; tree->SetBranchAddress("FatJet_deepTagMD_ZHbbvsQCD",&FatJet_deepTagMD_ZHbbvsQCD);
Float_t FatJet_deepTagMD_ZHccvsQCD[size_max]; tree->SetBranchAddress("FatJet_deepTagMD_ZHccvsQCD",&FatJet_deepTagMD_ZHccvsQCD);
Float_t FatJet_msoftdrop[size_max]; tree->SetBranchAddress("FatJet_msoftdrop",&FatJet_msoftdrop);

Float_t FatJet_deepTag_ZvsQCD[size_max]; tree->SetBranchAddress("FatJet_deepTag_ZvsQCD",&FatJet_deepTag_ZvsQCD);
UInt_t nFatJet = tree->SetBranchAddress("nFatJet",&nFatJet);
 Float_t Fatjet_pt[size_max]; tree->SetBranchAddress("FatJet_pt",&Fatjet_pt);
 Float_t Fatjet_mass[size_max]; tree->SetBranchAddress("FatJet_mass",&Fatjet_mass);
 Float_t Fatjet_phi[size_max]; tree->SetBranchAddress("FatJet_phi",&Fatjet_phi);
 Float_t Fatjet_eta[size_max]; tree->SetBranchAddress("FatJet_eta",&Fatjet_eta);
 TH1F *hist = new TH1F("hist","hist",30,0.,1.);


 TLorentzVector muon_vect, jet_vect, muon_vect_gen, b_vect,met_vect,el_vect_gen,muon_vect_boosted,jet_vect_boosted;
 float delta_R_true,delta_R_data, mass_inv_true, mass_inv_data, delta_R_true_corrected,mass_inv_true_corrected, angle, prod, deltaRTemp;
 double costhetap;
 
 TopUtilities Top1,top2; 


 float mu_pt_merged = 0.; 
 float mu_e_merged = 0.;
 float mu_phi_merged = 0.; 
 float mu_eta_merged = 0.; 
float mu_ch_merged = 0.;
//float mu_DB_merged = 0.;
float mu_Dxyerr_merged = 0.; 
float mu_Dz_merged = 0.; 
float mu_Dxy_merged = 0.;
float mu_MiniIso_merged = 0.;
float mu_Iso_merged = 0.;
float mu_IsGlobal_merged = 0.; 
float mu_IsTracker_merged = 0.;
float mu_NumberMatchedStations_merged = 0.;
float mu_Dxy_fract_merged=0.;
float pt_rel_merged=0.;
float mu_isHigh_merged=0.;
float mu_isTight_merged=0.;
int mu_high_truth_merged = 0;
int tau_high_truth_merged = 0;
int musize_merged = 0;
float nHadZ_merged =0;
float nLepTop_merged =0;

float nHadZ_merged_reco=0;
float nLepTop_merged_reco=0;

 float jet_pt_merged =0.; 
 float jet_e_merged = 0.;
 float jet_phi_merged = 0.; 
 float jet_eta_merged = 0.;
 int jet_high_truth_merged =0;
int jetsize_merged = 0;

 float top_pt_merged = 0.; 
 float top_e_merged = 0.;
 float top_phi_merged = 0.; 
 float top_eta_merged = 0.;
 float top_M_merged = 0.;
 int top_high_truth_merged = 0;
int topsize_merged = 0;


 float mu_boosted_pt_merged = 0.; 
 float mu_boosted_e_merged = 0.;
 float mu_boosted_phi_merged = 0.; 
 float mu_boosted_eta_merged = 0.;

  float jet_boosted_pt_merged = 0.; 
 float jet_boosted_e_merged = 0.;
 float jet_boosted_phi_merged = 0.; 
 float jet_boosted_eta_merged = 0.;

  float top_nu_pt_merged = 0.; 
 float top_nu_e_merged = 0.;
 float top_nu_phi_merged = 0.; 
 float top_nu_eta_merged = 0.;
 float top_nu_M_merged = 0.;

 double costheta_merged=0.;

 Float_t mu_pt_resolved = 0.; 
 float mu_e_resolved = 0.;
 float mu_phi_resolved = 0.; 
 float mu_eta_resolved = 0.; 
float mu_ch_resolved = 0.;
//float mu_DB_resolved = 0.;
float mu_Dxyerr_resolved = 0.; 
float mu_Dz_resolved = 0.; 
float mu_Dxy_resolved = 0.;
float mu_MiniIso_resolved = 0.;
float mu_Iso_resolved = 0.;
float mu_IsGlobal_resolved = 0.; 
float mu_IsTracker_resolved = 0.;
float mu_Dxy_fract_resolved=0.;
float pt_rel_resolved=0.;
float mu_isHigh_resolved=0.;
float mu_isTight_resolved=0.;
float mu_NumberMatchedStations_resolved = 0.;
int mu_high_truth_resolved = 0;
int tau_high_truth_resolved = 0;
int musize_resolved = 0;
float nHadZ_resolved =0;
float nLepTop_resolved =0;

float nHadZ_resolved_reco=0;
float nLepTop_resolved_reco=0;

 float jet_pt_resolved =0.; 
 float jet_e_resolved = 0.;
 float jet_phi_resolved = 0.; 
 float jet_eta_resolved = 0.;
 int jet_high_truth_resolved =0;
int jetsize_resolved = 0;

 float top_pt_resolved = 0.; 
 float top_e_resolved = 0.;
 float top_phi_resolved = 0.; 
 float top_eta_resolved = 0.;
 float top_M_resolved = 0.;
 int top_high_truth_resolved = 0;
int topsize_resolved = 0;

 float mu_boosted_pt_resolved = 0.; 
 float mu_boosted_e_resolved = 0.;
 float mu_boosted_phi_resolved = 0.; 
 float mu_boosted_eta_resolved = 0.;

  float jet_boosted_pt_resolved = 0.; 
 float jet_boosted_e_resolved = 0.;
 float jet_boosted_phi_resolved = 0.; 
 float jet_boosted_eta_resolved = 0.;

  float top_nu_pt_resolved = 0.; 
 float top_nu_e_resolved = 0.;
 float top_nu_phi_resolved = 0.; 
 float top_nu_eta_resolved = 0.;
 float top_nu_M_resolved = 0.;

 double costheta_resolved=0.;
bool pass_Htag, pass_Ztag, pass_Hadtag; 


TFile *f2= new TFile((output).c_str(),"RECREATE");
cout<<"begin"<<(output).c_str()<<endl;
TTree *is_merged= new TTree("is_merged","is_merged");
//branch 1 top merged
is_merged->Branch ( "Muon_Pt",&mu_pt_merged);
is_merged->Branch ( "Muon_Phi",&mu_phi_merged);
is_merged->Branch ( "Muon_Eta",&mu_eta_merged);
is_merged->Branch ( "Muon_E",&mu_e_merged);
is_merged->Branch ("Muon_Size",&musize_merged);
is_merged->Branch ("Muon_Charge",&mu_ch_merged);
///is_merged->Branch ("Muon_DB",&mu_DB_merged);
is_merged->Branch ("Muon_Dxyerr",&mu_Dxyerr_merged);
is_merged->Branch ("Muon_Dz",&mu_Dz_merged);
is_merged->Branch ("Muon_MiniIso",&mu_MiniIso_merged);
is_merged->Branch ("Muon_Iso",&mu_Iso_merged);
is_merged->Branch ("Muon_IsTrackerMuon",&mu_IsTracker_merged);
is_merged->Branch ("Muon_NumberMatchedStations",&mu_NumberMatchedStations_merged);
is_merged->Branch ("Muon_Dxy",&mu_Dxy_merged);
is_merged->Branch ("Muon_IsGlobalMuon",&mu_IsGlobal_merged);
is_merged->Branch ("Muon_High_Truth",&mu_high_truth_merged);
is_merged->Branch ("Muon_Dxy_fract",&mu_Dxy_fract_merged);
is_merged->Branch ("Muon_Pt_Rel",&pt_rel_merged);
is_merged->Branch ("Muon_isHigh",&mu_isHigh_merged);
is_merged->Branch ("Muon_isTight",&mu_isTight_merged);
is_merged->Branch ("Jet_High_Truth",&jet_high_truth_merged);
is_merged->Branch ("Top_High_Truth",&top_high_truth_merged);
is_merged->Branch ("Tau_High_Truth",&tau_high_truth_merged);

is_merged->Branch ( "Jet_Pt",&jet_pt_merged);
is_merged->Branch ( "Jet_Phi",&jet_phi_merged);
is_merged->Branch ( "Jet_Eta",&jet_eta_merged);
is_merged->Branch ( "Jet_E",&jet_e_merged);
is_merged->Branch ("Jet_Size",&jetsize_merged);

is_merged->Branch ( "top_Pt",&top_pt_merged);
is_merged->Branch ( "top_Phi",&top_phi_merged);
is_merged->Branch ( "top_Eta",&top_eta_merged);
is_merged->Branch ( "top_E",&top_e_merged);
is_merged->Branch ( "top_M",&top_M_merged);
is_merged->Branch ("top_Size",&topsize_merged);
is_merged->Branch ("Event_nHadZ",&nHadZ_merged);
is_merged->Branch ("Event_nLepTop",&nLepTop_merged);

is_merged->Branch ("Event_nHadZ_reco",&nHadZ_merged_reco);
is_merged->Branch ("Event_nLepTop_reco",&nLepTop_merged_reco);

is_merged->Branch ( "Muon_Boosted_Pt",&mu_boosted_pt_merged);
is_merged->Branch ( "Muon_Boosted_Phi",&mu_boosted_phi_merged);
is_merged->Branch ( "Muon_Boosted_Eta",&mu_boosted_eta_merged);
is_merged->Branch ( "Muon_Boosted_E",&mu_boosted_e_merged);

is_merged->Branch ( "Jet_Boosted_Pt",&jet_boosted_pt_merged);
is_merged->Branch ( "Jet_Boosted_Phi",&jet_boosted_phi_merged);
is_merged->Branch ( "Jet_Boosted_Eta",&jet_boosted_eta_merged);
is_merged->Branch ( "Jet_Boosted_E",&jet_boosted_e_merged);

is_merged->Branch ( "top_nu_Pt",&top_nu_pt_merged);
is_merged->Branch ( "top_nu_Phi",&top_nu_phi_merged);
is_merged->Branch ( "top_nu_Eta",&top_nu_eta_merged);
is_merged->Branch ( "top_nu_E",&top_nu_e_merged);
is_merged->Branch ( "top_nu_M",&top_nu_M_merged);

is_merged->Branch("Costheta",&costheta_merged);
is_merged->Branch("Pass_Htag",&pass_Htag);
is_merged->Branch("Pass_Hadtag",&pass_Hadtag);
is_merged->Branch("Pass_Ztag",&pass_Ztag);



TTree *is_resolved= new TTree("is_resolved","is_resolved");
//branch 1 top resolved
is_resolved->Branch ( "Muon_Pt",&mu_pt_resolved);
is_resolved->Branch ( "Muon_Phi",&mu_phi_resolved);
is_resolved->Branch ( "Muon_Eta",&mu_eta_resolved);
is_resolved->Branch ( "Muon_E",&mu_e_resolved);
is_resolved->Branch ("Muon_Size",&musize_resolved);
is_resolved->Branch ("Muon_Charge",&mu_ch_resolved);
//is_resolved->Branch ("Muon_DB",&mu_DB_resolved);
is_resolved->Branch ("Muon_Dxyerr",&mu_Dxyerr_resolved);
is_resolved->Branch ("Muon_Dz",&mu_Dz_resolved);
is_resolved->Branch ("Muon_MiniIso",&mu_MiniIso_resolved);
is_resolved->Branch ("Muon_Iso",&mu_Iso_resolved);
is_resolved->Branch ("Muon_IsTrackerMuon",&mu_IsTracker_resolved);
is_resolved->Branch ("Muon_NumberMatchedStations",&mu_NumberMatchedStations_resolved);
is_resolved->Branch ("Muon_Dxy",&mu_Dxy_resolved);
is_resolved->Branch ("Muon_IsGlobalMuon",&mu_IsGlobal_resolved);
is_resolved->Branch ("Muon_Dxy_fract",&mu_Dxy_fract_resolved);
is_resolved->Branch ("Muon_Pt_Rel",&pt_rel_resolved);
is_resolved->Branch ("Muon_isHigh",&mu_isHigh_resolved);
is_resolved->Branch ("Muon_isTight",&mu_isTight_resolved);
is_resolved->Branch ("Muon_High_Truth",&mu_high_truth_resolved);
is_resolved->Branch ("Jet_High_Truth",&jet_high_truth_resolved);
is_resolved->Branch ("Top_High_Truth",&top_high_truth_resolved);
is_resolved->Branch ("Tau_High_Truth",&tau_high_truth_resolved);

is_resolved->Branch ( "Jet_Pt",&jet_pt_resolved);
is_resolved->Branch ( "Jet_Phi",&jet_phi_resolved);
is_resolved->Branch ( "Jet_Eta",&jet_eta_resolved);
is_resolved->Branch ( "Jet_E",&jet_e_resolved);
is_resolved->Branch ("Jet_Size",&jetsize_resolved);

is_resolved->Branch ( "top_Pt",&top_pt_resolved);
is_resolved->Branch ( "top_Phi",&top_phi_resolved);
is_resolved->Branch ( "top_Eta",&top_eta_resolved);
is_resolved->Branch ( "top_E",&top_e_resolved);
is_resolved->Branch ( "top_M",&top_M_resolved);
is_resolved->Branch ("top_Size",&topsize_resolved);
is_resolved->Branch ("Event_nHadZ",&nHadZ_resolved);
is_resolved->Branch ("Event_nLepTop",&nLepTop_resolved);
is_resolved->Branch ("Event_nHadZ_reco",&nHadZ_resolved_reco);
is_resolved->Branch ("Event_nLepTop_reco",&nLepTop_resolved_reco);

is_resolved->Branch ( "Muon_Boosted_Pt",&mu_boosted_pt_resolved);
is_resolved->Branch ( "Muon_Boosted_Phi",&mu_boosted_phi_resolved);
is_resolved->Branch ( "Muon_Boosted_Eta",&mu_boosted_eta_resolved);
is_resolved->Branch ( "Muon_Boosted_E",&mu_boosted_e_resolved);

is_resolved->Branch ( "Jet_Boosted_Pt",&jet_boosted_pt_resolved);
is_resolved->Branch ( "Jet_Boosted_Phi",&jet_boosted_phi_resolved);
is_resolved->Branch ( "Jet_Boosted_Eta",&jet_boosted_eta_resolved);
is_resolved->Branch ( "Jet_Boosted_E",&jet_boosted_e_resolved);

is_resolved->Branch ( "top_nu_Pt",&top_nu_pt_resolved);
is_resolved->Branch ( "top_nu_Phi",&top_nu_phi_resolved);
is_resolved->Branch ( "top_nu_Eta",&top_nu_eta_resolved);
is_resolved->Branch ( "top_nu_E",&top_nu_e_resolved);
is_resolved->Branch ( "top_nu_M",&top_nu_M_resolved);

is_resolved->Branch("Costheta",&costheta_resolved);
is_resolved->Branch("Pass_Htag",&pass_Htag);
is_resolved->Branch("Pass_Hadtag",&pass_Hadtag);
is_resolved->Branch("Pass_Ztag",&pass_Ztag);
//is_resolved->Branch ("Muon_DB_fract",&mu_DB_fract_merged);
//is_resolved->Branch ("Muon_Pt_Rel",&pt_rel_merged);
 float Event_Number=0.;
float  Event_run=0.;


    int counter_merged=0, counter_resolved=0;
if(f!=0){
  nentries=f;
}
 for(int i=0; i<100000;i++){
  
  if(i%10000==1){
    cout<<i<<endl;

  }
  

  tree->GetEntry(i);

  mu_pt_resolved=0; 
  mu_e_resolved=0;
  mu_phi_resolved=0; 
  mu_eta_resolved=0;
  mu_ch_resolved=0; 
  
  mu_Dxyerr_resolved=0; 
  mu_Dz_resolved=0; 
  mu_Dxy_resolved=0; 
  mu_MiniIso_resolved=0;
  mu_Iso_resolved=0; 
  mu_IsGlobal_resolved=0; 
  mu_IsTracker_resolved=0; 
  mu_NumberMatchedStations_resolved=0; 
  mu_Dxy_fract_resolved=0.;
  pt_rel_resolved=0.;
  mu_isHigh_resolved=0.;
  mu_isTight_resolved=0.;
  mu_high_truth_resolved=0;
  jet_high_truth_resolved=0;
  top_high_truth_resolved=0;
  nHadZ_resolved=0;
  nLepTop_resolved=0;
  nHadZ_resolved_reco=0;
  nLepTop_resolved_reco=0;

  jet_pt_resolved=0; 
  jet_e_resolved=0;
  jet_phi_resolved=0; 
  jet_eta_resolved=0;

 top_pt_resolved=0; 
  top_e_resolved=0;
  top_phi_resolved=0; 
  top_eta_resolved=0;
  top_M_resolved=0;

   mu_boosted_pt_resolved = 0.; 
 mu_boosted_e_resolved = 0.;
 mu_boosted_phi_resolved = 0.; 
 mu_boosted_eta_resolved = 0.;

 tau_high_truth_resolved=0.;

  jet_boosted_pt_resolved = 0.; 
 jet_boosted_e_resolved = 0.;
 jet_boosted_phi_resolved = 0.; 
 jet_boosted_eta_resolved = 0.;

  top_nu_pt_resolved = 0.; 
 top_nu_e_resolved = 0.;
 top_nu_phi_resolved = 0.; 
 top_nu_eta_resolved = 0.;
 top_nu_M_resolved = 0.;

 costheta_resolved=0.;

   mu_pt_merged=0; 
  mu_e_merged=0;
  mu_phi_merged=0; 
  mu_eta_merged=0;
  mu_ch_merged=0; 

  mu_Dxy_fract_merged=0;
  mu_isHigh_merged=0.;
  mu_isTight_merged=0.;
  pt_rel_merged=0.; 
  mu_Dxyerr_merged=0; 
  mu_Dz_merged=0; 
  mu_Dxy_merged=0; 
  mu_MiniIso_merged=0;
   mu_Iso_merged=0;  
  mu_IsGlobal_merged=0; 
  mu_IsTracker_merged=0; 
  mu_NumberMatchedStations_merged=0; 
  mu_high_truth_merged=0;
  jet_high_truth_merged=0;
  top_high_truth_merged=0;
  tau_high_truth_merged=0.;
  nHadZ_merged=0;
  nLepTop_merged=0;
    nHadZ_merged_reco=0;
  nLepTop_merged_reco=0;

  jet_pt_merged=0; 
  jet_e_merged=0;
  jet_phi_merged=0; 
  jet_eta_merged=0;

 top_pt_merged=0; 
  top_e_merged=0;
  top_phi_merged=0; 
  top_eta_merged=0;
  top_M_merged=0;


   mu_boosted_pt_merged = 0.; 
 mu_boosted_e_merged = 0.;
 mu_boosted_phi_merged = 0.; 
 mu_boosted_eta_merged = 0.;

  jet_boosted_pt_merged = 0.; 
 jet_boosted_e_merged = 0.;
 jet_boosted_phi_merged = 0.; 
 jet_boosted_eta_merged = 0.;

  top_nu_pt_merged = 0.; 
 top_nu_e_merged = 0.;
 top_nu_phi_merged = 0.; 
 top_nu_eta_merged = 0.;
 top_nu_M_merged = 0.;

 costheta_merged=0.;
 
    float ntop_reco[100];
  float top_rec[100];
  int muon_high_truth[100];
  int top_high_truth[100];
   int jet_high_truth[100];
   int tau_high_truth[100];


int nTop_El=0;
int Z_tag=0;
int nZ_had=0;
int nTop_Lep=0;
int Top_tag_resolved=0;  
int Top_tag_merged=0;
met_vect.SetPtEtaPhiM(met_pt,0,met_phi,met_pt);

pass_Htag= false; 
pass_Ztag= false, 
pass_Hadtag= false;

for(int p=0;p<nFatJet;p++){
  if(FatJet_deepTag_H[p]>0.6 && FatJet_msoftdrop[p]>100 && FatJet_msoftdrop[p]<160) pass_Htag= true;
  if(FatJet_msoftdrop[p]>80 && FatJet_msoftdrop[p]<140 && FatJet_deepTag_ZvsQCD[p]>0.3) pass_Ztag= true;
}
  if(pass_Ztag==true  || pass_Htag==true ) pass_Hadtag=true;


  



for(int k=0;k<jetsize;k++){
  jet_high_truth[k]=0;
 


  if(jet_isDeep[k]>0.5 && met_pt>70 ){ //&& nTop_Lep>0


              jet_vect.SetPtEtaPhiM(jet_pt[k],jet_eta[k],jet_phi[k],jet_mass[k]);

              for(int u=0; u<gP_size;u++){

                  Int_t index1=0;
                 index1 = (gP_mom1I[u]);

                 
                 if(index1>-1){

                 if(abs(gP_flavor[u])==5 && (gP_flavor[u]*jet_fl[k])>0. && abs(gP_flavor[index1])==6 ){      // 

                 TLorentzVector jet_vect_gen;
                 jet_vect_gen.SetPtEtaPhiM(gP_Pt[u],gP_Eta[u],gP_Phi[u],gP_mass[u]);
                  float delta_R_true_corrected=sqrt(pow(gP_Eta[u]-jet_eta[k],2)+pow(jet_vect_gen.DeltaPhi(jet_vect),2));
                  if(delta_R_true_corrected<0.4){
                    jet_high_truth[k]=1;


                  }// MC match
                }// part giusta
              }
              }//loop generatore
        
          /*
              if(Jet_partonFlavour[k]==5 or Jet_partonFlavour[k]==-5){
                jet_high_truth[k]=1;
              }
*/
          for(int j=0; j<musize; j++){
            muon_high_truth[k]=0;
            tau_high_truth[k]=0;
            top_rec[k]=0;
            if(mu_pt[j]>10){
              muon_vect.SetPtEtaPhiM(mu_pt[j],mu_eta[j],mu_phi[j],mu_mass[j]);
              deltaRTemp=  sqrt(pow(jet_eta[k]-mu_eta[j],2)+pow(muon_vect.DeltaPhi(jet_vect),2));;

             
              if(deltaRTemp<2. && deltaRTemp>0.4  ){ //&& mu_isHigh[j]==1  && mu_Iso[j]<0.2
                top_high_truth[k]=0;
                top_rec[k]=1;
                TLorentzVector tot_vect_3=jet_vect+muon_vect;
                math::PtEtaPhiELorentzVector tot_vect= Top1.top4Momentum(muon_vect,jet_vect,met_vect.Px(),met_vect.Py());
                TLorentzVector tot_vect_1,tot_vect_2;
                tot_vect_1.SetPxPyPzE(-tot_vect.Px(),-tot_vect.Py(),-tot_vect.Pz(),tot_vect.E());

                top_nu_e_resolved=(tot_vect.E());
                top_nu_phi_resolved=(tot_vect.Phi());
                top_nu_pt_resolved=(tot_vect.Pt()); 
                top_nu_eta_resolved=(tot_vect.Eta()); 
                //top_reco_resolved=(top_rec[k]);
                top_nu_M_resolved=(sqrt(tot_vect.M2()));

                top_e_resolved=(tot_vect_3.E());
                top_phi_resolved=(tot_vect_3.Phi());
                top_pt_resolved=(tot_vect_3.Pt()); 
                top_eta_resolved=(tot_vect_3.Eta()); 

                top_M_resolved=(sqrt(tot_vect_3.M2()));
                //muon_index_resolved=(j);
                //jet_index_resolved=(k);
                jet_high_truth_resolved=jet_high_truth[k];
                tot_vect_2.SetPxPyPzE(tot_vect.Px(),tot_vect.Py(),tot_vect.Pz(),tot_vect.E());                
                costhetap= top2.costhetapol(muon_vect,jet_vect,tot_vect_2);

                mu_pt_resolved=mu_pt[j]; 
                mu_e_resolved=muon_vect.E();
                mu_phi_resolved=mu_phi[j]; 
                mu_eta_resolved=mu_eta[j];
                mu_ch_resolved=mu_ch[j];

                jet_pt_resolved=jet_pt[k]; 
                jet_e_resolved=jet_vect.E();
                jet_phi_resolved=jet_phi[k]; 
                jet_eta_resolved=jet_eta[k];

                 
                mu_Dxyerr_resolved=mu_Dxyerr[j]; 
                mu_Dz_resolved=mu_Dz[j]; 
                mu_Dxy_resolved=mu_Dxy[j]; 
                mu_MiniIso_resolved=mu_MiniIso[j];
                               mu_Iso_resolved=mu_Iso[j]; 
                mu_IsGlobal_resolved=mu_IsGlobal[j]; 
                mu_IsTracker_resolved=mu_IsTracker[j]; 
                mu_NumberMatchedStations_resolved=mu_NumberMatchedStations[j]; 
                mu_isHigh_resolved=mu_isHigh[j];
                mu_isTight_resolved=mu_isTight[j]; 
                mu_Dxy_fract_resolved=mu_Dxy[j]/(mu_Dxyerr[j]);
                pt_rel_resolved=((muon_vect.Vect()).Cross(jet_vect.Vect())).Mag()/((jet_vect.Vect()).Mag()); 
                for(int u=0; u<gP_size;u++){
                  int index1=gP_mom1I[u];
                 if(abs(gP_flavor[u])==13 && (gP_flavor[u]*mu_ch[j])<0. && abs(gP_flavor[index1])==24){       
                  muon_vect_gen.SetPtEtaPhiM(gP_Pt[u],gP_Eta[u],gP_Phi[u],gP_mass[u]);
                  float delta_R_true_corrected=sqrt(pow(gP_Eta[u]-mu_eta[j],2)+pow(muon_vect_gen.DeltaPhi(muon_vect),2));
                  if(delta_R_true_corrected<0.05){
                    muon_high_truth[k]=1;


                  }// MC match
                }// part giusta
               }//loop generatore

                for(int u=0; u<gP_size;u++){
                  int index1=gP_mom1I[u];
                 if(abs(gP_flavor[u])==15 && (gP_flavor[u]*mu_ch[j])<0. && abs(gP_flavor[index1])==24){       
                  muon_vect_gen.SetPtEtaPhiM(gP_Pt[u],gP_Eta[u],gP_Phi[u],gP_mass[u]);
                  float delta_R_true_corrected=sqrt(pow(gP_Eta[u]-mu_eta[j],2)+pow(muon_vect_gen.DeltaPhi(muon_vect),2));
                  if(delta_R_true_corrected<0.05){
                    tau_high_truth[k]=1;
                    

                  }// MC match
                }// part giusta
               }//loop generatore
               /*
                if( (int)mu_match[k] == 1){
                  muon_high_truth[k]=1;
                }
                else if( (int) mu_match[k]==15){
                  tau_high_truth[k]=1;
                }

                */
               if((int)mu_match[j] !=1 &&  muon_high_truth[k]==1){
                cout<<"fatal error"<<endl;
               }

              muon_vect.Boost(tot_vect_1.BoostVector());
              jet_vect.Boost(tot_vect_1.BoostVector());
 
              jet_boosted_pt_resolved=(jet_vect.Pt());
              jet_boosted_e_resolved=(jet_vect.E());
              jet_boosted_phi_resolved=(jet_vect.Phi());
              jet_boosted_eta_resolved=(jet_vect.Eta());

              mu_boosted_pt_resolved=(muon_vect.Pt());
              mu_boosted_e_resolved=(muon_vect.E());
              mu_boosted_phi_resolved=(muon_vect.Phi());
              mu_boosted_eta_resolved=(muon_vect.Eta());
              costheta_resolved=(costhetap);

              jet_vect.SetPtEtaPhiM(jet_pt[k],jet_eta[k],jet_phi[k],jet_mass[k]);
              muon_vect.SetPtEtaPhiM(mu_pt[j],mu_eta[j],mu_phi[j],mu_mass[j]);

              if((muon_high_truth[k]==1 || tau_high_truth[k]==1) && jet_high_truth[k]==1 && mu_ch[j]*jet_fl[k]>0){
                top_high_truth[k]=1;

              }
              mu_high_truth_resolved=(muon_high_truth[k]);
              tau_high_truth_resolved= tau_high_truth[k];
              top_high_truth_resolved=(top_high_truth[k]);
              counter_resolved=counter_resolved+1;
              is_resolved->Fill();
            }// selezione delta resolved

            else if(deltaRTemp<0.4){ // && mu_isHigh[j]==1 && abs( mu_DB[j]/(mu_DBerr[j]))<2.5 && (mu_pt[j]/(jet_pt[k]))>0.1
               top_high_truth[k]=0;
                top_rec[k]=-1;
                 TLorentzVector tot_vect_3;
                tot_vect_3=jet_vect;
                jet_vect=jet_vect-muon_vect;
               
                
               math::PtEtaPhiELorentzVector tot_vect= Top1.top4Momentum(muon_vect,jet_vect,met_vect.Px(),met_vect.Py());
                
                TLorentzVector tot_vect_1,tot_vect_2;
                tot_vect_1.SetPxPyPzE(-tot_vect.Px(),-tot_vect.Py(),-tot_vect.Pz(),tot_vect.E());
                top_nu_e_merged=(tot_vect.E());
                top_nu_phi_merged=(tot_vect.Phi());
                top_nu_pt_merged=(tot_vect.Pt()); 
                top_nu_eta_merged=(tot_vect.Eta()); 
                //top_reco_merged=(top_rec[k]);
                top_nu_M_merged=sqrt(tot_vect.M2());
                //muon_index_merged=(j);
                //jet_index_merged=(k);
                jet_high_truth_merged=jet_high_truth[k];
                top_e_merged=(tot_vect_3.E());
                top_phi_merged=(tot_vect_3.Phi());
                top_pt_merged=(tot_vect_3.Pt()); 
                top_eta_merged=(tot_vect_3.Eta()); 

                top_M_merged=sqrt(tot_vect_3.M2());
                tot_vect_2.SetPxPyPzE(tot_vect.Px(),tot_vect.Py(),tot_vect.Pz(),tot_vect.E());                
                costhetap= top2.costhetapol(muon_vect,jet_vect,tot_vect_2);   

                mu_pt_merged=mu_pt[j]; 
                mu_e_merged=muon_vect.E();
                mu_phi_merged=mu_phi[j]; 
                mu_eta_merged=mu_eta[j];
                mu_ch_merged=mu_ch[j];
                mu_Dxy_fract_merged=mu_Dxy[j]/(mu_Dxyerr[j]);
                pt_rel_merged=((muon_vect.Vect()).Cross(tot_vect_3.Vect())).Mag()/((tot_vect_3.Vect()).Mag());

                jet_pt_merged=jet_pt[k]; 
                jet_e_merged=jet_vect.E();
                jet_phi_merged=jet_phi[k]; 
                jet_eta_merged=jet_eta[k];

                //mu_D_merged=mu_DB[j]; 
                mu_Dxyerr_merged=mu_Dxyerr[j]; 
                mu_Dz_merged=mu_Dz[j]; 
                mu_Dxy_merged=mu_Dxy[j]; 
                mu_MiniIso_merged=mu_MiniIso[j];
                mu_Iso_merged=mu_Iso[j]; 
                mu_IsGlobal_merged=mu_IsGlobal[j]; 
                mu_IsTracker_merged=mu_IsTracker[j]; 
                mu_NumberMatchedStations_merged=mu_NumberMatchedStations[j];   
                mu_isHigh_merged=mu_isHigh[j];
                mu_isTight_merged=mu_isTight[j];      

                for(int u=0; u<gP_size;u++){
                  int index1=gP_mom1I[u];
                 if(abs(gP_flavor[u])==13 && (gP_flavor[u]*mu_ch[j])<0. && abs(gP_flavor[index1])==24){       
                  muon_vect_gen.SetPtEtaPhiM(gP_Pt[u],gP_Eta[u],gP_Phi[u],gP_mass[u]);
                  float delta_R_true_corrected=sqrt(pow(gP_Eta[u]-mu_eta[j],2)+pow(muon_vect_gen.DeltaPhi(muon_vect),2));
                  if(delta_R_true_corrected<0.1){
                    muon_high_truth[k]=1;


                  }// MC match
                }// part giusta
               }//loop generatore

                for(int u=0; u<gP_size;u++){
                  int index1=gP_mom1I[u];
                 if(abs(gP_flavor[u])==15 && (gP_flavor[u]*mu_ch[j])<0. && abs(gP_flavor[index1])==24){       
                  muon_vect_gen.SetPtEtaPhiM(gP_Pt[u],gP_Eta[u],gP_Phi[u],gP_mass[u]);
                  float delta_R_true_corrected=sqrt(pow(gP_Eta[u]-mu_eta[j],2)+pow(muon_vect_gen.DeltaPhi(muon_vect),2));
                  if(delta_R_true_corrected<0.05){
                    tau_high_truth[k]=1;
                    

                  }// MC match
                }// part giusta
               }//loop generatore
               
               /*
                if( (int)mu_match[k] == 1){
                  muon_high_truth[k]=1;
                }
                else if( (int) mu_match[k]==15){
                  tau_high_truth[k]=1;
                }

                */
              jet_vect=jet_vect+muon_vect;
              muon_vect.Boost(tot_vect_1.BoostVector());
              jet_vect.Boost(tot_vect_1.BoostVector());

              jet_boosted_pt_merged=(jet_vect.Pt());
              jet_boosted_e_merged=(jet_vect.E());
              jet_boosted_phi_merged=(jet_vect.Phi());
              jet_boosted_eta_merged=(jet_vect.Eta());

              mu_boosted_pt_merged=(muon_vect.Pt());
              mu_boosted_e_merged=(muon_vect.E());
              mu_boosted_phi_merged=(muon_vect.Phi());
              mu_boosted_eta_merged=(muon_vect.Eta());
              costheta_merged=(costhetap);


              jet_vect.SetPtEtaPhiM(jet_pt[k],jet_eta[k],jet_phi[k],jet_mass[k]);
              muon_vect.SetPtEtaPhiM(mu_pt[j],mu_eta[j],mu_phi[j],mu_mass[j]);

              if((muon_high_truth[k]==1 || tau_high_truth[k]==1) && jet_high_truth[k]==1 && mu_ch[j]*jet_fl[k]>0){
                top_high_truth[k]=1;

              }
              mu_high_truth_merged=(muon_high_truth[k]);
              top_high_truth_merged=(top_high_truth[k]);
              tau_high_truth_merged= tau_high_truth[k];
              counter_merged=counter_merged+1;
              is_merged->Fill();
             }// selezione delta merged
            }// muone giusto
            
          }// muon loop 
  }// b tagged
}//loop jet








}//end of loop entries

f1->Close();
is_merged->Write();

is_resolved->Write();






f2->Close();


cout<<"Done"<<(output).c_str()<<endl;
cout<<"counter_merged = "<<counter_merged<<endl;
cout<<"counter_resolved = "<<counter_resolved<<endl;

}// end of main 

void Tprime_final_conv_nano2(){

 Tprime_top1("Tprime_tHq_1000.root","Tprime_tHq_1000_ML.root",0);   

}


