#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <stdint.h>
#include <iomanip>
#include <ctime>
#include <map>
#include <math.h>

#include "TROOT.h"
#include "TFile.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TChain.h"
#include "TProfile.h"
#include "TMath.h"
#include "TString.h"
#include "TClass.h"
#include "TApplication.h"
#include "TLorentzVector.h"

#include "../interface/analysisUtils.h"
#include "../interface/setOutputTreeSynch.h"

using namespace std;

//*******MAIN*******************************************************************

int main (int argc, char** argv)
{ 
  std::string inputFolder = argv[1];
  std::string outputFile = argv[2];
  std::string leptonName = argv[3];
  std::string inputFile = argv[4];

  std::cout<<"file: "<<(inputFile).c_str()<<std::endl;
  TFile * fS = new TFile((inputFolder+leptonName+"/WZTree_"+inputFile+".root").c_str());
  TTree * inputTree = (TTree *)fS->Get("otree");

  int run;
  int event;
  int lumi;
  int njets;
  int nPV;
  int issignal;
  float pfMET;
  float pfMET_Phi;
  float l_pt;
  float l_eta;
  float l_phi;
  float l_e;
  float ungroomed_jet_pt;
  float ungroomed_jet_eta;
  float ungroomed_jet_phi;
  float ungroomed_jet_e;
  float jet_mass_pr;
  float jet_mass_so;
  float jet_tau2tau1;
  float v_pt;
  float v_eta;
  float v_phi;
  float v_mt;
  float mass_lvj_type0;
  int nBTagJet_medium;
  float jet2_pt;
  float jet2_btag;
  float jet3_pt;
  float jet3_btag;

  inputTree->SetBranchAddress("run", &run);
  inputTree->SetBranchAddress("event", &event);
  inputTree->SetBranchAddress("lumi", &lumi);
  inputTree->SetBranchAddress("njets", &njets);
  inputTree->SetBranchAddress("nPV", &nPV);
  inputTree->SetBranchAddress("issignal", &issignal);
  inputTree->SetBranchAddress("pfMET", &pfMET);
  inputTree->SetBranchAddress("pfMET_Phi", &pfMET_Phi);
  inputTree->SetBranchAddress("l_pt", &l_pt);
  inputTree->SetBranchAddress("l_eta", &l_eta);
  inputTree->SetBranchAddress("l_phi", &l_phi);
  inputTree->SetBranchAddress("l_e", &l_e);
  inputTree->SetBranchAddress("ungroomed_jet_pt", &ungroomed_jet_pt);
  inputTree->SetBranchAddress("ungroomed_jet_eta", &ungroomed_jet_eta);
  inputTree->SetBranchAddress("ungroomed_jet_phi", &ungroomed_jet_phi);
  inputTree->SetBranchAddress("ungroomed_jet_e", &ungroomed_jet_e);
  inputTree->SetBranchAddress("jet_mass_pr", &jet_mass_pr);
  inputTree->SetBranchAddress("jet_mass_so", &jet_mass_so);
  inputTree->SetBranchAddress("jet_tau2tau1", &jet_tau2tau1);
  inputTree->SetBranchAddress("v_pt", &v_pt);
  inputTree->SetBranchAddress("v_eta", &v_eta);
  inputTree->SetBranchAddress("v_phi", &v_phi);
  inputTree->SetBranchAddress("v_mt", &v_mt);
  inputTree->SetBranchAddress("mass_lvj_type2", &mass_lvj_type0);
  inputTree->SetBranchAddress("nBTagJet_medium", &nBTagJet_medium);
  inputTree->SetBranchAddress("jet2_pt", &jet2_pt);
  inputTree->SetBranchAddress("jet2_btag", &jet2_btag);
  inputTree->SetBranchAddress("jet3_pt", &jet3_pt);
  inputTree->SetBranchAddress("jet3_btag", &jet3_btag);

  //---------output tree----------------
  TFile* outROOT = TFile::Open((std::string("output/output_synch_")+leptonName+std::string("/")+std::string("WZTree_")+outputFile+std::string(".root")).c_str(),"recreate");
  outROOT->cd();
  TTree* outTree = new TTree("otree", "otree");
  outTree->SetDirectory(0);
  setOutputTreeSynch *WZTree = new setOutputTreeSynch(outTree);
  
  //---------start loop on events------------
  for (Long64_t jentry=0; jentry<inputTree->GetEntries();jentry++) {

    inputTree->GetEntry(jentry);

    WZTree->initializeVariables(); //initialize all variables

    if(jentry % 1000 == 0)    
      cout << "read entry: " << jentry << endl;

    WZTree->issignal = issignal;

    //save event variables
    WZTree->run   = run;
    WZTree->event = event;
    WZTree->lumi = lumi;
    WZTree->njets = njets;
    WZTree->nPV  = nPV;
    
    WZTree->l_pt  = l_pt;
    WZTree->l_eta = l_eta;
    WZTree->l_phi = l_phi;	

    WZTree->pfMET   = pfMET;
    WZTree->pfMETPhi = pfMET_Phi;
    
    WZTree->W_pt = v_pt;
    WZTree->W_eta = v_eta;
    WZTree->W_phi = v_phi;

    WZTree->jet_pt  = ungroomed_jet_pt;
    WZTree->jet_eta = ungroomed_jet_eta;
    WZTree->jet_phi = ungroomed_jet_phi;
    WZTree->jet_mass_pruned   = jet_mass_pr;
    WZTree->jet_mass_softdrop   = jet_mass_so;
    WZTree->jet_tau2tau1   = jet_tau2tau1;

    WZTree->nbtag=nBTagJet_medium;
    WZTree->m_lvj = mass_lvj_type0;

    WZTree->jet2_pt = jet2_pt;
    WZTree->jet2_btag = jet2_btag;
    WZTree->jet3_pt = jet3_pt;
    WZTree->jet3_btag = jet3_btag;

    //fill the tree
    if(strcmp(leptonName.c_str(),"mu")==0 && WZTree->issignal==1 && WZTree->W_pt>200 && WZTree->pfMET>40 && WZTree->l_pt>53 && WZTree->jet_pt>200 && WZTree->nbtag <1 && ((WZTree->jet_mass_pruned > 40 && WZTree->jet_mass_pruned<65) || (WZTree->jet_mass_pruned > 135 && WZTree->jet_mass_pruned<150))) {
      //WZTree->jet_mass_pruned > 40 && WZTree->jet_mass_pruned < 130) {//&& WZTree->jet_tau2tau1<0.5) {
      outTree->Fill();
    }
    if(strcmp(leptonName.c_str(),"el")==0 && WZTree->issignal==1 && WZTree->W_pt>200 && WZTree->pfMET>80 && WZTree->l_pt>120 && WZTree->jet_pt>200 && WZTree->nbtag <1 && ((WZTree->jet_mass_pruned > 40 && WZTree->jet_mass_pruned<65) || (WZTree->jet_mass_pruned > 135 && WZTree->jet_mass_pruned<150))) {
      //WZTree->jet_mass_pruned > 40 && WZTree->jet_mass_pruned < 130) {//&& WZTree->jet_tau2tau1<0.5) {
      outTree->Fill();
    }
  }

  //--------close everything-------------
  outTree->Write();
  outROOT->Close();

  return(0);
}
