#include "../interface/setOutputTree.h"

int run;
int event;
int nJets;
int nVtx;
float met;
float met_px;
float met_py;
float met_pz_type0;
float met_pz_type2;
float leptonPt;
float leptonEta;
float leptonPhi;
float leptonE;
float AK8jetPt[10];
float AK8jetEta[10];
float AK8jetPhi[10];
float AK8jetE[10];
float AK8jetPrunedMass[10];
float AK8jetTrimmedMass[10];
float AK8jetFilteredMass[10];
float AK8jetTau1[10];
float AK8jetTau2[10];
float AK8jetTau3[10];
float jetPt[10];
float jetEta[10];
float jetPhi[10];
float jetE[10];
float jet_bDiscr[10];
int genBosonPdgId[10];
float genBosonPt[10];
float genBosonEta[10];
float genBosonPhi[10];
float genBosonE[10];
int genLeptonPdgId[10];
float genLeptonPt[10];
float genLeptonEta[10];
float genLeptonPhi[10];
float genLeptonE[10];
float genNuPt[10];
float genNuEta[10];
float genNuPhi[10];
float genNuE[10];
float deltaR_lak8jet;
float deltaphi_METak8jet;
float deltaphi_Vak8jet;
float W_pt;
float W_eta;
float W_phi;
float W_E;
float W_mt;
float boosted_lvj_m_type0;
float boosted_lvj_m_type2;

// List of branches
TBranch *b_run;
TBranch *b_event;
TBranch *b_nJets;
TBranch *b_nVtx;
TBranch *b_met;
TBranch *b_met_px;
TBranch *b_met_py;
TBranch *b_met_pz_type0;
TBranch *b_met_pz_type2;
TBranch *b_leptonPt;
TBranch *b_leptonEta;
TBranch *b_leptonPhi;
TBranch *b_leptonE;
TBranch *b_AK8jetPt;
TBranch *b_AK8jetEta;
TBranch *b_AK8jetPhi;
TBranch *b_AK8jetE;
TBranch *b_AK8jetPrunedMass;
TBranch *b_AK8jetTrimmedMass;
TBranch *b_AK8jetFilteredMass;
TBranch *b_AK8jetTau1;
TBranch *b_AK8jetTau2;
TBranch *b_AK8jetTau3;
TBranch *b_jetPt;
TBranch *b_jetEta;
TBranch *b_jetPhi;
TBranch *b_jetE;
TBranch *b_jet_bDiscr;
TBranch *b_genBosonPdgId;
TBranch *b_genBosonPt;
TBranch *b_genBosonEta;
TBranch *b_genBosonPhi;
TBranch *b_genBosonE;
TBranch *b_genLeptonPdgId;
TBranch *b_genLeptonPt;
TBranch *b_genLeptonEta;
TBranch *b_genLeptonPhi;
TBranch *b_genLeptonE;
TBranch *b_genNuPt;
TBranch *b_genNuEta;
TBranch *b_genNuPhi;
TBranch *b_genNuE;
TBranch *b_deltaR_lak8jet;
TBranch *b_deltaphi_METak8jet;
TBranch *b_deltaphi_Vak8jet;
TBranch *b_W_pt;
TBranch *b_W_eta;
TBranch *b_W_phi;
TBranch *b_W_E;
TBranch *b_W_mt;
TBranch *b_boosted_lvj_m_type0;
TBranch *b_boosted_lvj_m_type2;

void init()
{
  run=-999;
  event=-999;
  nJets=-999;
  nVtx=-999;
  met=-999;
  met_px=-999;
  met_py=-999;
  met_pz_type0=-999;
  met_pz_type2=-999;
  leptonPt=-999;
  leptonEta=-999;
  leptonPhi=-999;
  leptonE=-999;
  deltaR_lak8jet=-999;
  deltaphi_METak8jet=-999;
  deltaphi_Vak8jet=-999;
  W_pt=-999;
  W_eta=-999;
  W_phi=-999;
  W_E=-999;
  W_mt=-999;
  boosted_lvj_m_type0=-999;
  boosted_lvj_m_type2=-999;

 for (int i=0; i<10; i++) {
   AK8jetPt[i]=-999;
   AK8jetEta[i]=-999;
   AK8jetPhi[i]=-999;
   AK8jetE[i]=-999;
   AK8jetPrunedMass[i]=-999;
   AK8jetTrimmedMass[i]=-999;
   AK8jetFilteredMass[i]=-999;
   AK8jetTau1[i]=-999;
   AK8jetTau2[i]=-999;
   AK8jetTau3[i]=-999;
   jetPt[i]=-999;
   jetEta[i]=-999;
   jetPhi[i]=-999;
   jetE[i]=-999;
   jet_bDiscr[i]=-999;
   genBosonPdgId[i]=-999;
   genBosonPt[i]=-999;
   genBosonEta[i]=-999;
   genBosonPhi[i]=-999;
   genBosonE[i]=-999;
   genLeptonPdgId[i]=-999;
   genLeptonPt[i]=-999;
   genLeptonEta[i]=-999;
   genLeptonPhi[i]=-999;
   genLeptonE[i]=-999;
   genNuPt[i]=-999;
   genNuEta[i]=-999;
   genNuPhi[i]=-999;
   genNuE[i]=-999;
 }
}

void SetOutTree(TTree* outTree)
{
  outTree->Branch("run",&event,"run/I");
  outTree->Branch("event",&event,"event/I");
  outTree->Branch("nJets",&nJets,"nJets/I");
  outTree->Branch("nVtx",&nVtx,"nVtx/I");
  outTree->Branch("met",&met,"met/F");
  outTree->Branch("met_px",&met_px,"met_px/F");
  outTree->Branch("met_py",&met_py,"met_py/F");
  outTree->Branch("met_pz_type0",&met_pz_type0,"met_pz_type0/F");
  outTree->Branch("met_pz_type2",&met_pz_type2,"met_pz_type2/F");
  outTree->Branch("leptonPt",&leptonPt,"leptonPt/F");
  outTree->Branch("leptonEta",&leptonEta,"leptonEta/F");
  outTree->Branch("leptonPhi",&leptonPhi,"leptonPhi/F");
  outTree->Branch("leptonE",&leptonE,"leptonE/F");
  outTree->Branch("AK8jetPt",&AK8jetPt,"AK8jetPt[10]/F");
  outTree->Branch("AK8jetEta",&AK8jetEta,"AK8jetEta[10]/F");
  outTree->Branch("AK8jetPhi",&AK8jetPhi,"AK8jetPhi[10]/F");
  outTree->Branch("AK8jetE",&AK8jetE,"AK8jetE[10]/F");
  outTree->Branch("AK8jetPrunedMass",&AK8jetPrunedMass,"AK8jetPrunedMass[10]/F");
  outTree->Branch("AK8jetTrimmedMass",&AK8jetTrimmedMass,"AK8jetTrimmedMass[10]/F");
  outTree->Branch("AK8jetFilteredMass",&AK8jetFilteredMass,"AK8jetFilteredMass[10]/F");
  outTree->Branch("AK8jetTau1",&AK8jetTau1,"AK8jetTau1[10]/F");
  outTree->Branch("AK8jetTau2",&AK8jetTau2,"AK8jetTau2[10]/F");
  outTree->Branch("AK8jetTau3",&AK8jetTau3,"AK8jetTau3[10]/F");
  outTree->Branch("jetPt",&jetPt,"jetPt[10]/F");
  outTree->Branch("jetEta",&jetEta,"jetEta[10]/F");
  outTree->Branch("jetPhi",&jetPhi,"jetPhi[10]/F");
  outTree->Branch("jetE",&jetE,"jetE[10]/F");
  outTree->Branch("jet_bDiscr",&jet_bDiscr,"jet_bDiscr[10]/F");
  outTree->Branch("genBosonPdgId",&genBosonPdgId,"genBosonPdgId[10]/I");
  outTree->Branch("genBosonPt",&genBosonPt,"genBosonPt[10]/F");
  outTree->Branch("genBosonEta",&genBosonEta,"genBosonEta[10]/F");
  outTree->Branch("genBosonPhi",&genBosonPhi,"genBosonPhi[10]/F");
  outTree->Branch("genBosonE",&genBosonE,"genBosonE[10]/F");
  outTree->Branch("genLeptonPdgId",&genLeptonPdgId,"genLeptonPdgId[10]/I");
  outTree->Branch("genLeptonPt",&genLeptonPt,"genLeptonPt[10]/F");
  outTree->Branch("genLeptonEta",&genLeptonEta,"genLeptonEta[10]/F");
  outTree->Branch("genLeptonPhi",&genLeptonPhi,"genLeptonPhi[10]/F");
  outTree->Branch("genLeptonE",&genLeptonE,"genLeptonE[10]/F");
  outTree->Branch("genNuPt",&genNuPt,"genNuPt[10]/F");
  outTree->Branch("genNuEta",&genNuEta,"genNuEta[10]/F");
  outTree->Branch("genNuPhi",&genNuPhi,"genNuPhi[10]/F");
  outTree->Branch("genNuE",&genNuE,"genNuE[10]/F");
  outTree->Branch("deltaR_lak8jet",&deltaR_lak8jet,"deltaR_lak8jet/F");
  outTree->Branch("deltaphi_METak8jet",&deltaphi_METak8jet,"deltaphi_METak8jet/F");
  outTree->Branch("deltaphi_Vak8jet",&deltaphi_Vak8jet,"deltaphi_Vak8jet/F");
  outTree->Branch("W_pt",&W_pt,"W_pt/F");
  outTree->Branch("W_eta",&W_eta,"W_eta/F");
  outTree->Branch("W_phi",&W_phi,"W_phi/F");
  outTree->Branch("W_E",&W_E,"W_E/F");
  outTree->Branch("W_mt",&W_mt,"W_mt/F");
  outTree->Branch("boosted_lvj_m_type0",&boosted_lvj_m_type0,"boosted_lvj_m_type0/F");
  outTree->Branch("boosted_lvj_m_type2",&boosted_lvj_m_type2,"boosted_lvj_m_type2/F");
}

void InitRecoTree(TTree* nt)
{
  nt->SetBranchAddress("run", &run, &b_run);
  nt->SetBranchAddress("event", &event, &b_event);
  nt->SetBranchAddress("met", &met, &b_met);
  nt->SetBranchAddress("met_px", &met_px, &b_met_px);
  nt->SetBranchAddress("met_py", &met_py, &b_met_py);
  nt->SetBranchAddress("met_pz_type0", &met_pz_type0, &b_met_pz_type0);
  nt->SetBranchAddress("met_pz_type2", &met_pz_type2, &b_met_pz_type2);
  nt->SetBranchAddress("leptonPt",&leptonPt,&b_leptonPt);
  nt->SetBranchAddress("leptonEta",&leptonEta,&b_leptonEta);
  nt->SetBranchAddress("leptonPhi",&leptonPhi,&b_leptonPhi);
  nt->SetBranchAddress("leptonE",&leptonE,&b_leptonE);
  nt->SetBranchAddress("AK8jetPt",&AK8jetPt,&b_AK8jetPt);
  nt->SetBranchAddress("AK8jetEta",&AK8jetEta,&b_AK8jetEta);
  nt->SetBranchAddress("AK8jetPhi",&AK8jetPhi,&b_AK8jetPhi);
  nt->SetBranchAddress("AK8jetE",&AK8jetE,&b_AK8jetE);
  nt->SetBranchAddress("AK8jetPrunedMass",&AK8jetPrunedMass,&b_AK8jetPrunedMass);
  nt->SetBranchAddress("AK8jetTrimmedMass",&AK8jetTrimmedMass,&b_AK8jetTrimmedMass);
  nt->SetBranchAddress("AK8jetFilteredMass",&AK8jetFilteredMass,&b_AK8jetFilteredMass);
  nt->SetBranchAddress("AK8jetTau1",&AK8jetTau1,&b_AK8jetTau1);
  nt->SetBranchAddress("AK8jetTau2",&AK8jetTau2,&b_AK8jetTau2);
  nt->SetBranchAddress("AK8jetTau3",&AK8jetTau3,&b_AK8jetTau3);
  nt->SetBranchAddress("jetPt",&jetPt,&b_jetPt);
  nt->SetBranchAddress("jetEta",&jetEta,&b_jetEta);
  nt->SetBranchAddress("jetPhi",&jetPhi,&b_jetPhi);
  nt->SetBranchAddress("jetE",&jetE,&b_jetE);
  nt->SetBranchAddress("jet_bDiscr",&jet_bDiscr,&b_jet_bDiscr);
  nt->SetBranchAddress("genBosonPdgId",&genBosonPdgId,&b_genBosonPdgId);
  nt->SetBranchAddress("genBosonPt",&genBosonPt,&b_genBosonPt);
  nt->SetBranchAddress("genBosonEta",&genBosonEta,&b_genBosonEta);
  nt->SetBranchAddress("genBosonPhi",&genBosonPhi,&b_genBosonPhi);
  nt->SetBranchAddress("genBosonE",&genBosonE,&b_genBosonE);
  nt->SetBranchAddress("genLeptonPdgId",&genLeptonPdgId,&b_genLeptonPdgId);
  nt->SetBranchAddress("genLeptonPt",&genLeptonPt,&b_genLeptonPt);
  nt->SetBranchAddress("genLeptonEta",&genLeptonEta,&b_genLeptonEta);
  nt->SetBranchAddress("genLeptonPhi",&genLeptonPhi,&b_genLeptonPhi);
  nt->SetBranchAddress("genLeptonE",&genLeptonE,&b_genLeptonE);
  nt->SetBranchAddress("genNuPt",&genNuPt,&b_genNuPt);
  nt->SetBranchAddress("genNuEta",&genNuEta,&b_genNuEta);
  nt->SetBranchAddress("genNuPhi",&genNuPhi,&b_genNuPhi);
  nt->SetBranchAddress("genNuE",&genNuE,&b_genNuE);
  nt->SetBranchAddress("deltaR_lak8jet",&deltaR_lak8jet,&b_deltaR_lak8jet);
  nt->SetBranchAddress("deltaphi_METak8jet",&deltaphi_METak8jet,&b_deltaphi_METak8jet);
  nt->SetBranchAddress("deltaphi_Vak8jet",&deltaphi_Vak8jet,&b_deltaphi_Vak8jet);
  nt->SetBranchAddress("W_pt",&W_pt,&b_W_pt);
  nt->SetBranchAddress("W_eta",&W_eta,&b_W_eta);
  nt->SetBranchAddress("W_phi",&W_phi,&b_W_phi);
  nt->SetBranchAddress("W_E",&W_E,&b_W_E);
  nt->SetBranchAddress("W_mt",&W_mt,&b_W_mt);
  nt->SetBranchAddress("boosted_lvj_m_type0",&boosted_lvj_m_type0,&b_boosted_lvj_m_type0);  
  nt->SetBranchAddress("boosted_lvj_m_type2",&boosted_lvj_m_type2,&b_boosted_lvj_m_type2);  
}
