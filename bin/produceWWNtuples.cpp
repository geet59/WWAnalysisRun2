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
#include "TF1.h"
#include <TClonesArray.h>           

#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TVertex.hh"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "BaconAna/DataFormats/interface/TAddJet.hh"
#include "BaconAna/Utils/interface/TTrigger.hh"
#include "BaconAna/DataFormats/interface/TLHEWeight.hh"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"


#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
#include "../BtagUnc.hh"

#include "../interface/setOutputTree.h"
#include "../interface/METzCalculator.h"
#include "../interface/METzCalculator_Run2.h"
#include "../interface/analysisUtils.h"
#include "../interface/readJSONFile.h"
#include "../interface/Utils.hh"

using namespace std;

//*******MAIN*******************************************************************

int main (int argc, char** argv)
{ 

  int t0 = time(NULL);

  std::string inputFolder = argv[1];
  std::string outputFile = argv[2];
  int isMC = atoi(argv[3]);
  std::string cluster = argv[4];
  std::string inputTreeName = argv[5];
  std::string inputFile = argv[6];
  std::string xSecWeight = argv[7];
  std::string TotalNumberOfEntries = argv[8];
  std::string TotalNumberOfNegativeEntries = argv[9];
  int applyTrigger = atoi(argv[10]);
  std::string jsonFileName = argv[11];
  int isLocal = atoi(argv[12]);
  int VBFSel  = atoi(argv[13]);

  std::string leptonName;

  if ( VBFSel==1)	cout<<"==> VBF selection method : Select two highest pT jets"<<endl;
  else if ( VBFSel==2)	cout<<"==> VBF selection method : Select pair with highest mjj..."<<endl;
  else if ( VBFSel==3)	cout<<"==> VBF selection method : Select pair with highest DeltaEta..."<<endl;
 // else if ( VBFSel==4)	cout<<"==> VBF selection method : Select two highest pT jets..."<<endl;
  else {	cout<<"\n\nERROR:	Enter valid vbf selection criteria....\n\n"<<endl;
    exit(0);  
  }

  std::string iHLTFile="${CMSSW_BASE}/src/BaconAna/DataFormats/data/HLTFile_25ns";
  const std::string cmssw_base = getenv("CMSSW_BASE");
  std::string cmssw_base_env = "${CMSSW_BASE}";
  size_t start_pos = iHLTFile.find(cmssw_base_env);
  if(start_pos != std::string::npos) {
    iHLTFile.replace(start_pos, cmssw_base_env.length(), cmssw_base);
  }

  const baconhep::TTrigger triggerMenu(iHLTFile);  
  std::cout<<"Apply trigger: "<<applyTrigger<<std::endl;


  TLorentzVector LEP1, LEP2, LEP3;// SJ1_PuppiAK8, SJ2_PuppiAK8, SJ1, SJ2;
  TLorentzVector NU0,NU1,NU2,NU3;//NU0_puppi,NU1_puppi,NU2_puppi;
  TLorentzVector NU0_jes_up, NU0_jes_dn;
  TLorentzVector NU0_jer_up, NU0_jer_dn;
  TLorentzVector  AK4;
  TLorentzVector JET_jes_up, JET_jes_dn;// JET_PuppiAK8_jes_up, JET_PuppiAK8_jes_dn;
  TLorentzVector AK4_JET1,AK4_JET2;
  TLorentzVector AK4_JET1_jes_up, AK4_JET1_jes_dn;
  TLorentzVector AK4_JET2_jes_up, AK4_JET2_jes_dn;
  TLorentzVector VBF1,VBF2,TOT;
  TLorentzVector VBF1_jes_up, VBF1_jes_dn, VBF2_jes_up, VBF2_jes_dn;
  TLorentzVector ELE,MU;


  std::vector<TLorentzVector> tightMuon;
  std::vector<TLorentzVector> looseMuon;
  std::vector<TLorentzVector> tightEle;
  std::vector<TLorentzVector> looseEle;

  int ok=0, total=0;
  baconhep::TEventInfo *info  	= new baconhep::TEventInfo();
  baconhep::TGenEventInfo *gen	= new baconhep::TGenEventInfo();
  TClonesArray *genPartArr 	= new TClonesArray("baconhep::TGenParticle");
  TClonesArray *electronArr	= new TClonesArray("baconhep::TElectron");
  TClonesArray *vertexArr	= new TClonesArray("baconhep::TVertex");
  TClonesArray *jetArr		= new TClonesArray("baconhep::TJet");
  TClonesArray *lheWgtArr	= new TClonesArray("baconhep::TLHEWeight");


/*  char command1[3000];
  char command2[3000];

  if ( cluster == "lxplus")
    sprintf(command1, "eos find -f %s  | awk '!/log|fail/ {print $1}' | awk 'NF {print \"root://eoscms.cern.ch/\"$1}' > listTemp_%s.txt", (inputFolder).c_str(), outputFile.c_str());	// NF in awk command skips the blank line
  else
sprintf(command1,"eos root://cmseos.fnal.gov find %s | grep root | awk '{print \"root://cmseos.fnal.gov/\"$1}' > listTemp_%s.txt",(inputFolder).c_str(),  outputFile.c_str());
 
   // sprintf(command1,"xrdfs root://cmseos.fnal.gov ls %s | awk '{print \"root://cmseos.fnal.gov/\"$1}' > listTemp_%s.txt",(inputFolder).c_str(),  outputFile.c_str());
  std::cout<<command1<<std::endl;
  sprintf(command2,"sed -i '/failed$/d' listTemp_%s.txt",outputFile.c_str());
  system(command1);
  system(command2);
  char list1[2000];
  sprintf (list1, "listTemp_%s.txt", outputFile.c_str());
  ifstream rootList (list1);
  char command3[300];
  sprintf(command3, "rm listTemp_%s.txt", outputFile.c_str());
  system(command3);
*/
 

char list1[2000];
   sprintf (list1, "%s", inputFile.c_str());
ifstream rootList (list1);

 int fileCounter=0;

  vector<TString>  sampleName; 

  while (!rootList.eof())
  {
    char iRun_tW[700];
    rootList >> iRun_tW;
    if(!rootList.good())break;
    sampleName.push_back(iRun_tW);
    fileCounter++;
  }

  TFile *infile=0;
  TTree *eventTree=0;
  int cutEff[21]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  TFile* pileupFileMC = TFile::Open("puWeights_80x_37ifb.root");
  TH1D* puWeights = (TH1D*)pileupFileMC->Get("puWeights");
  TH1D* puWeightsUp = (TH1D*)pileupFileMC->Get("puWeightsUp");
  TH1D* puWeightsDown = (TH1D*)pileupFileMC->Get("puWeightsDown");
  puWeights->SetBins(75,0,75);
  puWeightsUp->SetBins(75,0,75);
  puWeightsDown->SetBins(75,0,75);
  TFile* IDIsoEle = TFile::Open("egammaEffi_EGM2D_TightCutBasedIDSF.root","READ");
  TH1F *hIDIsoEle = (TH1F*)IDIsoEle->Get("EGamma_SF2D");

  TFile* GSFCorrEle = TFile::Open("egammaEffi_SF2D_GSF_tracking.root","READ");
  TH1F *hGSFCorrEle = (TH1F*)GSFCorrEle->Get("EGamma_SF2D");

  TFile* TriggerEle = TFile::Open("ElectronTrigger_SF.root","READ");
  TH1F* hTriggerEle = (TH1F*)TriggerEle->Get("HLT_Ele27");

  TFile* IDMuA = TFile::Open("MuonID_RunBCDEF_23SepReReco_19p72fb.root","READ");
  TH1F *hIDMuA = (TH1F*)IDMuA->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio");

  TFile* IDMuB = TFile::Open("MuonID_RunGH_23SepReReco_16p146fb.root","READ");
  TH1F *hIDMuB = (TH1F*)IDMuB->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio");

  TFile* IsoMuA = TFile::Open("MuonIso_RunBCDEF_23SepReReco_19p72fb.root","READ");
  TH1F *hIsoMuA = (TH1F*)IsoMuA->Get("TightISO_TightID_pt_eta/abseta_pt_ratio");

  TFile* IsoMuB = TFile::Open("MuonIso_RunGH_23SepReReco_16p146fb.root","READ");
  TH1F *hIsoMuB = (TH1F*)IsoMuB->Get("TightISO_TightID_pt_eta/abseta_pt_ratio");

  TFile* TriggerMuA = TFile::Open("MuonTrigger_RunBCDEF_23SepReReco_19p72fb.root","READ");
  TH1F* hTriggerMuA = (TH1F*)TriggerMuA->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio");

  TFile* TriggerMuB = TFile::Open("MuonTrigger_RunGH_23SepReReco_16p146fb.root","READ");
  TH1F* hTriggerMuB = (TH1F*)TriggerMuB->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio");

  TFile* file = TFile::Open( "puppiCorr.root","READ");
  TFile* outROOT = TFile::Open((outputFile+(".root")).c_str(),"recreate");
  outROOT->cd();
  TTree* outTree = new TTree("otree", "otree");
  setOutputTree* WZTree = new setOutputTree(outTree);

  int nEvents=0;	
  Long64_t jentry2=0;
  int count_genEvents=0;

  int nInputFiles = sampleName.size();

  if (isLocal==1) nInputFiles = 1;
  cout<<"==> Total number of input files : "<<nInputFiles<<endl;

  TH1D *MCpu = new TH1D("MCpu","",75,0,75);
  TH1D *MCpu_up = new TH1D("MCpu_up","",75,0,75);
  TH1D *MCpu_down = new TH1D("MCpu_down","",75,0,75);

  Long64_t TotalNumberOfEvents = 0, nNegEvents = 0; 

  BTagCalibration calib("csvv2", "CSVv2_Moriond17_B_H.csv");
  BTagCalibrationReader bTagReader(BTagEntry::OP_LOOSE,  // working point: can be OP_LOOSE, OP_MEDIUM, OP_TIGHT 
      "central",             // label for the central value (see the scale factor file)
      {"up","down"});        // vector of labels for systematics
  bTagReader.load(calib, BTagEntry::FLAV_B, "comb");      // use the "comb" measurements for b-jets
  bTagReader.load(calib, BTagEntry::FLAV_C, "comb");      // use the "comb" measurements for c-jets
  bTagReader.load(calib, BTagEntry::FLAV_UDSG, "incl");   // use the "incl" measurements for light jets


  for(int i=0;i<nInputFiles;i++)
  {//loop10 begins
    infile = TFile::Open(sampleName[i]);
    eventTree = (TTree*)infile->Get("Events");

    TotalNumberOfEvents+=eventTree->GetEntries();
    if(isMC)
    { 
      eventTree->SetBranchAddress("Info", &info);    TBranch *infoBr = eventTree->GetBranch("Info");
      TBranch *genBr=0;
      eventTree->SetBranchAddress("GenEvtInfo", &gen); genBr = eventTree->GetBranch("GenEvtInfo");
      for (Long64_t jentry=0; jentry<eventTree->GetEntries();jentry++,jentry2++)
      {
	//eventTree->GetEntry(jentry);
	genBr->GetEntry(jentry);
	infoBr->GetEntry(jentry);	    
	MCpu->Fill(info->nPUmean);
	MCpu_up->Fill(info->nPUmeanp);
	MCpu_down->Fill(info->nPUmeanm);
	if (jentry2%50000 == 0) std::cout << "\t File no. " << i << "; Neg Event Count; read entry: " << jentry2 <<"/"<<TotalNumberOfEvents<<std:: endl;
	if (gen->weight<0)	nNegEvents++;
      }
    }
    delete infile;
    infile=0, eventTree=0;
  }  ///loop10 ends


  cout<<"==> Total number of events : "<<TotalNumberOfEvents<<endl;
  cout<<"==> Total number of negative events : "<<nNegEvents<<endl;

  //  float weight = std::atof(xSecWeight.c_str())/((TotalNumberOfEvents) - 2*std::atof(TotalNumberOfNegativeEntries.c_str()));
  float weight = std::atof(xSecWeight.c_str())/((TotalNumberOfEvents) - 2*std::atof(TotalNumberOfNegativeEntries.c_str()));
  cout<<"Weight of cross-sec/events = "<<weight<<endl;
  int totalEntries=0;


  JetCorrectorParameters paramAK4puppi("Summer16_23Sep2016V4_MC_JEC/Summer16_23Sep2016V4_MC_Uncertainty_AK4PFPuppi.txt");
  JetCorrectorParameters paramAK4chs("Summer16_23Sep2016V4_MC_JEC/Summer16_23Sep2016V4_MC_Uncertainty_AK4PFchs.txt");
  JetCorrectorParameters paramAK8chs("Summer16_23Sep2016V4_MC_JEC/Summer16_23Sep2016V4_MC_Uncertainty_AK8PFchs.txt");
  JetCorrectorParameters paramAK8puppi("Summer16_23Sep2016V4_MC_JEC/Summer16_23Sep2016V4_MC_Uncertainty_AK8PFPuppi.txt");
  JetCorrectionUncertainty *fJetUnc_AK4chs = new JetCorrectionUncertainty(paramAK4chs);
  //---------start loop on events------------
  std::cout << "---------start loop on events------------" << std::endl;
  jentry2=0;
  for(int i=0;i<nInputFiles;i++)
  {
    cout<<"\n\n=====	Processing File Number : "<<i<<"/"<<nInputFiles<<"\n\t"<<sampleName[i]<<"\n-------"<<endl;

    infile = TFile::Open(sampleName[i]);
    eventTree = (TTree*)infile->Get("Events");

    totalEntries+=eventTree->GetEntries();

    nEvents=eventTree->GetEntries();

    cout<<"\t==> Entries = "<<nEvents<<endl;



    eventTree->SetBranchAddress("Info", &info);    TBranch *infoBr = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("Electron", &electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");
    eventTree->SetBranchAddress("PV",   &vertexArr); TBranch *vertexBr = eventTree->GetBranch("PV");
    eventTree->SetBranchAddress("AK4CHS",   &jetArr); TBranch *jetBr = eventTree->GetBranch("AK4CHS");    
    TBranch *genBr=0, *genPartBr=0, *lhePartBr=0;
    if(isMC)
    {  //loop11 begins 
      eventTree->SetBranchAddress("GenEvtInfo", &gen); genBr = eventTree->GetBranch("GenEvtInfo");
      eventTree->SetBranchAddress("GenParticle",&genPartArr); genPartBr = eventTree->GetBranch("GenParticle");
      if(eventTree->GetListOfBranches()->FindObject("LHEWeight"))
      {
	eventTree->SetBranchAddress("LHEWeight",&lheWgtArr); lhePartBr = eventTree->GetBranch("LHEWeight");	       }
    } //loop11 ends
    for (Long64_t jentry=0; jentry<eventTree->GetEntries();jentry++,jentry2++)
    {
      infoBr->GetEntry(jentry);	    
//if	(jentry2>100) exit(0);


      int GenPassCut = 0;

      tightMuon.clear();
      tightEle.clear();
      looseMuon.clear();
      looseEle.clear();
      if (jentry2%10000 == 0) std::cout << "\tread entry: " << jentry2 <<"/"<<TotalNumberOfEvents<<std:: endl;
      WZTree->initializeVariables(); //initialize all variables

      WZTree->run   = info->runNum;
      WZTree->event = info->evtNum;
      WZTree->lumi  = info->lumiSec;

      if (isMC==1)
      {//loop14 begins
	lheWgtArr->Clear();
	if(lhePartBr)
	{
	  lhePartBr->GetEntry(jentry);
	}
	genPartArr->Clear();
	genBr->GetEntry(jentry);
	genPartBr->GetEntry(jentry);
	TLorentzVector hadW, lepW, VBFJ1, VBFJ2, VBFJ, temp;
	TLorentzVector genLep, genNeutrino, genWquarks, genVBFquarks;
	std::vector<TLorentzVector> v_GEN_hadW, v_GEN_lepW, v_GEN_VBFJ1, v_GEN_VBFJ2, v_GEN_VBFJ, v_GEN_temp;
	std::vector<TLorentzVector> vZ_genLep,vW_genLep, v_genNeutrino, v_genWquarks, v_genVBFquarks;

	v_GEN_hadW.clear();	v_GEN_lepW.clear();	v_GEN_VBFJ1.clear();	v_GEN_VBFJ2.clear();
	v_GEN_VBFJ.clear();	v_GEN_temp.clear();	vZ_genLep.clear();	v_genNeutrino.clear();
	v_genWquarks.clear(); 	vW_genLep.clear();//	v_genVBFquarks.clear();
	for (int i = 0; i<genPartArr->GetEntries();i++)
	{  //loop12 begins
	  const baconhep::TGenParticle* genloop = (baconhep::TGenParticle*) ((*genPartArr)[i]);
	  Int_t parentPdg=dynamic_cast<baconhep::TGenParticle *>(genPartArr->At(genloop->parent>-1 ? genloop->parent : 0))->pdgId;
	  if( (abs(genloop->pdgId) == 11 || abs(genloop->pdgId) == 13 ) && abs(parentPdg) == 23)
	  {
	    genLep.SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
	    vZ_genLep.push_back(genLep);
	  }
	  if( (abs(genloop->pdgId) == 11 || abs(genloop->pdgId) == 13 ) && abs(parentPdg) == 24)
	  {
	    genLep.SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
	    vW_genLep.push_back(genLep);
	  }
	  if( (abs(genloop->pdgId) == 12 || abs(genloop->pdgId) == 14 ) && abs(parentPdg) == 24)
	  {
	    genNeutrino.SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
	    v_genNeutrino.push_back(genNeutrino);
	  }

	  if( (abs(genloop->pdgId) == 1 || abs(genloop->pdgId) == 3 || abs(genloop->pdgId) == 5 || abs(genloop->pdgId) == 2 || abs(genloop->pdgId) == 4 || abs(genloop->pdgId) == 6) && genloop->status == 23 && abs(parentPdg) != 24)
	  {
	    genVBFquarks.SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
	    v_genVBFquarks.push_back(genVBFquarks);
	  }
	} //loop12 ends 
	if (vZ_genLep.size()==2 && vW_genLep.size()==1 && v_genNeutrino.size()==1 && v_genVBFquarks.size()==2)
	{ //loop13 begins
	  WZTree->isGen           = 1;
	  WZTree->lep1_pt_gen      = (vZ_genLep[0]).Pt();
	  WZTree->lep1_eta_gen     = (vZ_genLep[0]).Eta();
	  WZTree->lep1_m_gen       = (vZ_genLep[0]).M();
	  WZTree->lep2_pt_gen      = (vZ_genLep[1]).Pt();
	  WZTree->lep2_eta_gen     = (vZ_genLep[1]).Eta();
	  WZTree->lep2_m_gen     = (vZ_genLep[1]).M();
	  WZTree->lep3_pt_gen      = (vW_genLep[2]).Pt();
	  WZTree->lep3_eta_gen     = (vW_genLep[2]).Eta();
	  WZTree->lep3_m_gen     = (vW_genLep[2]).M();
	  WZTree->nu3_pz_gen	  = (v_genNeutrino[2]).Pz();
	  WZTree->nu3_pt_gen	  = (v_genNeutrino[2]).Pt();
	  WZTree->nu3_eta_gen	  = (v_genNeutrino[2]).Eta();


	  WZTree->dilep_m_gen    = (vZ_genLep[0]+vZ_genLep[1]).M();

	  WZTree->trilep_m_gen    = (vZ_genLep[0]+vZ_genLep[1]+vW_genLep[2]).M();    


	  WZTree->lepZ_pt_gen     = (vZ_genLep[0]+vZ_genLep[1]).Pt();
	  WZTree->lepZ_eta_gen    = (vZ_genLep[0]+vZ_genLep[1]).Eta();
	  WZTree->lepZ_phi_gen    = (vZ_genLep[0]+vZ_genLep[1]).Phi();
	  WZTree->lepZ_e_gen      = (vZ_genLep[0]+vZ_genLep[1]).E();
	  WZTree->lepZ_m_gen      = (vZ_genLep[0]+vZ_genLep[1]).M();
	  WZTree->lepW_pt_gen     = (vW_genLep[2]+v_genNeutrino[2]).Pt();
	  WZTree->lepW_eta_gen    = (vW_genLep[2]+v_genNeutrino[2]).Eta();
	  WZTree->lepW_phi_gen    = (vW_genLep[2]+v_genNeutrino[2]).Phi();
	  WZTree->lepW_e_gen      = (vW_genLep[2]+v_genNeutrino[2]).E();
	  WZTree->lepW_m_gen      = (vW_genLep[2]+v_genNeutrino[2]).M();

	  WZTree->WZ_mass_gen	= (vZ_genLep[0] + vZ_genLep[1] + vW_genLep[2] + v_genNeutrino[2]).M();
	  WZTree->WZ_mT_gen	= (vZ_genLep[0] + vZ_genLep[1] + vW_genLep[2] + v_genNeutrino[2]).Mt();
	  WZTree->WZ_pT_gen	= (vZ_genLep[0] + vZ_genLep[1] + vW_genLep[2] + v_genNeutrino[2]).Pt();
	  WZTree->WZ_eta_gen    = (vZ_genLep[0] + vZ_genLep[1] + vW_genLep[2] + v_genNeutrino[2]).Eta();


	  WZTree->AK4_1_pt_gen	= v_genVBFquarks[0].Pt();
	  WZTree->AK4_1_eta_gen	= v_genVBFquarks[0].Eta();
	  WZTree->AK4_1_phi_gen	= v_genVBFquarks[0].Phi();
	  WZTree->AK4_1_e_gen	= v_genVBFquarks[0].E();
	  WZTree->AK4_1_mass_gen	= v_genVBFquarks[0].M();

	  WZTree->AK4_2_pt_gen	= v_genVBFquarks[1].Pt();
	  WZTree->AK4_2_eta_gen	= v_genVBFquarks[1].Eta();
	  WZTree->AK4_2_phi_gen	= v_genVBFquarks[1].Phi();
	  WZTree->AK4_2_e_gen	= v_genVBFquarks[1].E();
	  WZTree->AK4_2_mass_gen	= v_genVBFquarks[1].M();

	  WZTree->AK4_jj_DeltaEta_gen = abs(v_genVBFquarks[0].Eta() - v_genVBFquarks[1].Eta());
	  WZTree->AK4_jj_mass_gen = (v_genVBFquarks[0] + v_genVBFquarks[1]).M();

	  WZTree->Zeppen1=abs(WZTree->WZ_eta_gen -(WZTree->AK4_1_eta_gen + WZTree->AK4_2_eta_gen)/2) ;




	  count_genEvents++;
	}
	if ( WZTree->lep1_pt_gen > 25 && WZTree->lep2_pt_gen >15 && WZTree->lep3_pt_gen >20 && ((abs(WZTree->lep1_eta_gen)) <2.5) && ((abs(WZTree->lep2_eta_gen)) <2.5) && ((abs(WZTree->lep3_eta_gen)) <2.5) && WZTree->dilep_m_gen >4 && WZTree->trilep_m_gen >100  && ((abs(WZTree->AK4_jj_DeltaEta_gen))>2.5) && WZTree->AK4_1_pt_gen >50 && WZTree->AK4_2_pt_gen>50 && ((abs(WZTree->AK4_1_eta_gen)) <4.7) && ((abs(WZTree->AK4_2_eta_gen)) <4.7) && WZTree->AK4_jj_mass_gen > 500   && ((abs((WZTree->dilep_m_gen) -91.2))<15) && WZTree->Zeppen1<2.5 ) 

	{
	  GenPassCut = 1;
	}
	for (int i = 0; i<lheWgtArr->GetEntries();i++)     // Note that i is starting from 446.
	{
	  const baconhep::TLHEWeight *lhe = (baconhep::TLHEWeight*)((*lheWgtArr)[i]);
	  WZTree->LHEWeight[i] = lhe->weight;
	}
      }//loop14 ends

      WZTree->issignal = 0;
      WZTree->wSampleWeight = weight; //xsec/TotalNumberOfEntries
      WZTree->pu_Weight = 1.; //temporary value
      WZTree->pu_Weight_up = 1.; //temporary value
      WZTree->pu_Weight_down = 1.; //temporary value
      WZTree->top1_NNLO_Weight = 1.;
      WZTree->top2_NNLO_Weight = 1.;
      WZTree->id_eff_Weight = 1.;
      WZTree->id_eff_Weight2 = 1.;
      WZTree->trig_eff_Weight = 1.;
      WZTree->trig_eff_Weight2 = 1.;
WZTree->id_eff_Weight3 = 1.;
      WZTree->trig_eff_Weight3 = 1.;

      if (gen->weight>0)
	WZTree->genWeight=1.;
      else if (gen->weight<0) {  //loop15 begins
	WZTree->genWeight=-1.;
	//nNegEvents++;
      }  //loop15 ends
      cutEff[0]++;

      if (isMC==1)
      {   //loop16 begins
	if (GenPassCut == 1)   
	  cutEff[1]++;
      } //loop16ends


      vertexArr->Clear();
      vertexBr->GetEntry(jentry);
      WZTree->nPV = vertexArr->GetEntries();
      //PILE-UP WEIGHT
      if (isMC==1) { //loop17 begins
	if(int(info->nPUmean)<75){
	  WZTree->pu_Weight = puWeights->GetBinContent(info->nPUmean); //our pu recipe
	  WZTree->pu_Weight_up = puWeightsUp->GetBinContent(info->nPUmean); //our pu recipe
	  WZTree->pu_Weight_down = puWeightsDown->GetBinContent(info->nPUmean); //our pu recipe
	}
	else{ //should not happen as we have a weight for all simulated n_pu multiplicities!
	  std::cout<<"Warning! n_pu too big"<<std::endl;
	  // throw logic_error("n_pu too big");
	  WZTree->pu_Weight = 0.;
	  WZTree->pu_Weight_up = 0.;
	  WZTree->pu_Weight_down = 0.;
	} 
      }//loop17 ends

      if(applyTrigger==1)
	if(! triggerMenu.pass("HLT_Ele27_WPTight_Gsf_v*",info->triggerBits)) continue;

      /////////////////THE SELECTED LEPTON
      int nTightEle=0, nLooseEle=0;
      double pt_cut = 15;//pvrevious=25
      double leadelept_cut = 15;
      electronArr->Clear();
      electronBr->GetEntry(jentry);
      const baconhep::TElectron *leadele = NULL;
      const baconhep::TElectron *subele = NULL;
      const baconhep::TElectron *subsubele = NULL;

      //double leadeleE=-999, subeleE=-999, subsubeleE=-999;
      double iso = 1.5;
      for (int i=0; i<electronArr->GetEntries(); i++) {    //loop1 on electron begins
	const baconhep::TElectron *ele = (baconhep::TElectron*)((*electronArr)[i]);
	//	cout<<jentry	<<"		" <<electronArr->GetEntries() <<"		"<<	i	<<"		"<<ele->pt<<endl;
	if (ele->pt<=pt_cut) continue;
	if (fabs(ele->eta)>=2.5) continue;
	if(!passEleLooseSel(ele,info->rhoIso)) continue;
	nLooseEle++;
	if(!passEleTightSel(ele,info->rhoIso)) continue;
	ELE.SetPtEtaPhiM(ele->pt,ele->eta,ele->phi,0.0005109989461);
	tightEle.push_back(ELE);
	nTightEle++;
	iso = ele->chHadIso + TMath::Max( 0.0,(ele->gammaIso + ele->neuHadIso - info->rhoIso*eleEffArea(ele->eta)) );
	if(!leadele || ele->pt>leadele->pt)
	{
	  if(!(ele->pt>leadelept_cut)) continue;
	  subele = leadele;
	  leadele = ele;
	  //leadeleE = ELE.E();
	  WZTree->l_iso1 = iso/ele->pt;
	}
	else if (!subele || ele->pt > subele->pt)
	{
	  subsubele =subele;
	  subele = ele;
	  //subeleE = ELE.E();
	  WZTree->l_iso2 = iso/ele->pt;
	}
	else if (!subsubele || ele->pt > subsubele->pt)
	{
	  subsubele = ele;
	  //subsubeleE = ELE.E();
	  WZTree->l_iso3 = iso/ele->pt;
	}

      } //loop1 on electron ends

      if (nTightEle != 3) continue;
      if (nLooseEle>3) continue;
      cutEff[2]++;
      //cout << jentry2 << "\tTLV: tightEle = " << tightEle.size() << "\t" << tightEle[0].Pt() << "\t" << tightEle[1].Pt() << "\t" << tightEle[2].Pt() << endl;
      if ((leadele->q>0 && subele->q>0 && subsubele>0) || (leadele->q<0 && subele->q<0 && subsubele<0))
          continue;
      cutEff[3]++;
      float charge1 = 0.0;
      float charge2 = 0.0;
      int tightEle_ZL1_index = -1;
      int tightEle_ZL2_index = -1;
      float deltaM = 900.0;
      for(int i=0; i<nTightEle-1; i++) {
          for(int j=i+1; j<nTightEle; j++) {
	      if (i==0) charge1 = leadele->q;
	      if (j==0) charge2 = leadele->q;
	      if (i==1) charge1 = subele->q;
	      if (j==1) charge2 = subele->q;
	      if (i==2) charge1 = subsubele->q;
	      if (j==2) charge2 = subsubele->q;
	      if (charge1*charge2 > 0) continue;
	      if (fabs((tightEle[i]+tightEle[j]).M()-91.1876) < deltaM) {
	          deltaM = fabs((tightEle[i]+tightEle[j]).M()-91.1876); 
		  tightEle_ZL1_index = i;
		  tightEle_ZL2_index = j;
              }
          }
      }
      cutEff[4]++;
      //Get index for w-boson
      int tightEle_WL_index = -1;
      for(int i=0; i<nTightEle; i++) {
          if (i==tightEle_ZL1_index || i==tightEle_ZL2_index) continue;
	  tightEle_WL_index = i;
      }
      //cout << jentry2 << "\t" << tightEle_ZL1_index << "\t" << tightEle_ZL2_index << "\t" << (tightEle[tightEle_ZL1_index]+tightEle[tightEle_ZL2_index]).M() << endl;

      // continue if m_ll is less then 4 GeV
//      if ((tightEle[tightEle_ZL1_index]+tightEle[tightEle_ZL2_index]).M() < 4.0) continue;
      cutEff[5]++;
      // contiue if Z-mass is outside mass window of 15 GeV
      if (fabs((tightEle[tightEle_ZL1_index]+tightEle[tightEle_ZL2_index]).M() - 91.1876) > 15.) continue;
      cutEff[6]++;
      // continue if m_lll is less then 100 GeV
    //  if ((tightEle[tightEle_ZL1_index]+tightEle[tightEle_ZL2_index]+tightEle[tightEle_WL_index]).M() < 100.) continue;
      //cutEff[7]++;

      //cout << "=====> " << jentry2 << "\t" << tightEle_ZL1_index << "\t" << tightEle_ZL2_index << "\t" << (tightEle[tightEle_ZL1_index]+tightEle[tightEle_ZL2_index]).M() << endl;

      WZTree->l_pt1 = tightEle[tightEle_ZL1_index].Pt();
      WZTree->l_eta1 = tightEle[tightEle_ZL1_index].Eta();
      WZTree->l_phi1 = tightEle[tightEle_ZL1_index].Phi();
      WZTree->l_e1 = tightEle[tightEle_ZL1_index].M();
      //WZTree->l_charge1 = chargel1;

      WZTree->l_pt2  = tightEle[tightEle_ZL2_index].Pt();
      WZTree->l_eta2 = tightEle[tightEle_ZL2_index].Eta();
      WZTree->l_phi2 = tightEle[tightEle_ZL2_index].Phi();
      WZTree->l_e2 = tightEle[tightEle_ZL2_index].M();
      //WZTree->l_charge2 = chargel2;

      WZTree->l_pt3  = tightEle[tightEle_WL_index].Pt();
      WZTree->l_eta3 = tightEle[tightEle_WL_index].Eta();
      WZTree->l_phi3 = tightEle[tightEle_WL_index].Phi();
      WZTree->l_e3 = tightEle[tightEle_WL_index].M();
      WZTree->dilep_m = (tightEle[tightEle_ZL1_index]+tightEle[tightEle_ZL2_index]).M();
      WZTree->trilep_m = (tightEle[tightEle_ZL1_index]+tightEle[tightEle_ZL2_index]+tightEle[tightEle_WL_index]).M();
//WZTree->l_charge3 = chargel3;

      // //preselection on met
      if (info->pfMETC < 30) continue;   //Et(miss)>40GeV
      cutEff[7]++;

      jetArr->Clear();
      jetBr->GetEntry(jentry);
      std::vector<int> indexGoodVBFJets;
      for ( int i=0; i<jetArr->GetEntries(); i++) 
      {const baconhep::TJet *jet = (baconhep::TJet*)((*jetArr)[i]);
	      bool isCleaned = true;
 
	      if (jet->pt<30) continue;
	      for ( std::size_t j=0; j<tightEle.size(); j++) {
		      if (deltaR(tightEle.at(j).Eta(), tightEle.at(j).Phi(),
        jet->eta,   jet->phi) < 0.4) {
			      isCleaned = false;
		      }
	      }
	     for ( std::size_t j=0; j<tightMuon.size(); j++) {
		      if (deltaR(tightMuon.at(j).Eta(), tightMuon.at(j).Phi(),
		      jet->eta,   jet->phi) < 0.4) {
			      isCleaned = false;
		      }
	     }
	     if (isCleaned==false) continue;
	      if (jet->pt<30) continue;
	      if (!passJetLooseSel(jet)) continue; 
	      if (fabs(jet->eta)>=4.7) continue;
	      indexGoodVBFJets.push_back(i);
	      WZTree->njets++;
	      AK4.SetPtEtaPhiM(jet->pt,jet->eta,jet->phi,jet->mass);
      }
      if (indexGoodVBFJets.size() <=2) continue;
      cutEff[8]++;

      int nVBF1=-1, nVBF2=-1; 
      double jetselectid[2]={-999, -999};
      double jetselectpt[2]={0, 0};
      for (std::size_t i=0; i<indexGoodVBFJets.size(); i++){
	      const baconhep::TJet *jet = (baconhep::TJet*)((*jetArr)[indexGoodVBFJets.at(i)]);
	      if(jet->pt<30) continue;
	      if(jet->pt>jetselectpt[1]){
		      jetselectpt[0]=jetselectpt[1];
		      jetselectid[0]=jetselectid[1];
		      jetselectpt[1]=jet->pt;
		      jetselectid[1]=i;
	      }
	      else if(jet->pt>jetselectpt[0])
	      {
		      jetselectpt[0]=jet->pt;
		      jetselectid[0]=i;
	      }
      }
      if (jetselectid[0]<0)continue;
      if (jetselectid[1]<0)continue;

	nVBF1 = indexGoodVBFJets.at(jetselectid[1]); //save position of the 1st vbf jet
	nVBF2 = indexGoodVBFJets.at(jetselectid[0]); //save position of the 2nd vbf jet

	const baconhep::TJet *jet1 = (baconhep::TJet*)((*jetArr)[nVBF1]);
	const baconhep::TJet *jet2 = (baconhep::TJet*)((*jetArr)[nVBF2]);
	VBF1.SetPtEtaPhiM(jet1->pt,jet1->eta,jet1->phi,jet1->mass);
	VBF2.SetPtEtaPhiM(jet2->pt,jet2->eta,jet2->phi,jet2->mass);
	TOT = VBF1 + VBF2;
	WZTree->vbf_maxpt_j1_pt = jet1->pt;
	WZTree->vbf_maxpt_j1_eta = jet1->eta;
	WZTree->vbf_maxpt_j1_phi = jet1->phi;
	WZTree->vbf_maxpt_j1_e = VBF1.E();
	WZTree->vbf_maxpt_j1_mass = VBF1.M();

	WZTree->vbf_maxpt_j2_pt = jet2->pt;
	WZTree->vbf_maxpt_j2_eta = jet2->eta;
	WZTree->vbf_maxpt_j2_phi = jet2->phi;
	WZTree->vbf_maxpt_j2_e = VBF2.E();
	WZTree->vbf_maxpt_j2_mass = VBF2.M();

	if (TOT.M()<500) continue;
	cutEff[9]++;
	WZTree->vbf_maxpt_jj_Deta = abs(VBF1.Eta() - VBF2.Eta());
	if (abs(VBF1.Eta() - VBF2.Eta()) <2.5) continue;
	cutEff[10]++;
	indexGoodVBFJets.clear();
	cutEff[11]++;
	WZTree->ZeppenfeldW1 =(((LEP1+LEP2+LEP3).Eta()) - ((VBF1.Eta() + VBF2.Eta())/2.0));

if(nTightEle>0){
      WZTree->type=1;
      leptonName = "el";	// Added this part for neutrino pz calculation in case there is w-boson.
}	if (strcmp(leptonName.c_str(),"el")==0 && isMC==1) {//loop 3 begins
	//  apply ID, ISO SF's
	WZTree->id_eff_Weight = GetSFs_Lepton(WZTree->l_pt1, WZTree->l_eta1, hIDIsoEle);	// Get Scale factor corresponding to the pt and eta.
	
	WZTree->id_eff_Weight = WZTree->id_eff_Weight*GetSFs_Lepton(WZTree->l_pt1, WZTree->l_eta1, hGSFCorrEle);
	WZTree->trig_eff_Weight = 1.0/GetSFs_Lepton(WZTree->l_pt1, WZTree->l_eta1, hTriggerEle);
	
	WZTree->id_eff_Weight2 = GetSFs_Lepton(WZTree->l_pt2, WZTree->l_eta2, hIDIsoEle);	// Get Scale factor corresponding to the pt and eta.
	
  						// apply GSF/RECO SF's for electrons 
	WZTree->id_eff_Weight2 = WZTree->id_eff_Weight2*GetSFs_Lepton(WZTree->l_pt2, WZTree->l_eta2, hGSFCorrEle);
	WZTree->trig_eff_Weight2 = 1.0/GetSFs_Lepton(WZTree->l_pt2, WZTree->l_eta2, hTriggerEle);
	
	WZTree->id_eff_Weight3 = GetSFs_Lepton(WZTree->l_pt3, WZTree->l_eta3, hIDIsoEle);       // Get Scale factor corresponding to the pt and eta.
	
	// apply GSF/RECO SF's for electrons
	WZTree->id_eff_Weight3 = WZTree->id_eff_Weight3*GetSFs_Lepton(WZTree->l_pt3, WZTree->l_eta3, hGSFCorrEle);
	WZTree->trig_eff_Weight3 = 1.0/GetSFs_Lepton(WZTree->l_pt3, WZTree->l_eta3, hTriggerEle);
	}
	cutEff[12]++;  												//loop 3 ends
      outTree->Fill();
      //cout<<"DEBUG: 2:" << endl;
      //cout<<"DEBUG: 3:" << endl;
    }//loop on entries end
    delete infile;
    infile=0, eventTree=0;
    /////////////////FILL THE TREE
  }//loop on file ends
  //delete puWeight;	delete puWeight_up;	delete puWeight_down;
  delete MCpu;	delete MCpu_up;	delete MCpu_down;
  delete puWeightsDown;	delete puWeightsUp;	delete puWeights;


  pileupFileMC->Close();
  file->Close();
  std::cout << "---------end loop on events------------" << std::endl;
  std::cout << std::endl;
  std::cout << "GEN events = " << count_genEvents << std::endl;



  std::cout << "----------------------" << std::endl;
  std::cout << " SUMMARY" << std::endl;
  std::cout << "----------------------" << std::endl;
  std::cout << std::endl;
  std::cout<<"MC matching: "<<(float)ok/(float)total<<std::endl;
  std::cout<<"negative events: "<<nNegEvents<<std::endl;
  std::cout << std::endl;
  std::cout<<"(0) all events:        "<<cutEff[0]<<"\t:\t"<<((float)cutEff[0]*100.0)/(float)cutEff[0]<<std::endl
	  <<"(1) Gen Events:        "<<cutEff[1]<<"\t:\t"<<((float)cutEff[1]*100.0)/(float)cutEff[0]<<std::endl
	  <<"(2) Exactly 3 electron:  "<<cutEff[2]<<"\t:\t"<<((float)cutEff[2]*100.0)/(float)cutEff[0]<<std::endl
	  <<"(3) charges of all electron not same:      "<<cutEff[3]<<"\t:\t"<<((float)cutEff[3]*100.0)/(float)cutEff[2]<<std::endl
	  <<"(4) two electron from Z:               "<<cutEff[4]<<"\t:\t"<<((float)cutEff[4]*100.0)/(float)cutEff[3]<<std::endl
	  <<"(5) 3rd electron from W:               "<<cutEff[5]<<"\t:\t"<<((float)cutEff[5]*100.0)/(float)cutEff[4]<<std::endl
	  <<"(6) Dilepton-Zmass:  "<<cutEff[6]<<"\t:\t"<<((float)cutEff[6]*100.0)/(float)cutEff[5]<<std::endl
	  <<"(7) MET:               "<<cutEff[7]<<"\t:\t"<<((float)cutEff[7]*100.0)/(float)cutEff[6]<<std::endl
	  <<"(8) >=2 good VBF jets: "<<cutEff[8]<<"\t:\t"<<((float)cutEff[8]*100.0)/(float)cutEff[7]<<std::endl
	  <<"(9) two Highest pt VBF jets with mjj>500     : "<<cutEff[9]<<"\t:\t"<<((float)cutEff[9]*100.0)/(float)cutEff[8]<<std::endl
	  <<"(10) Events with VBFjets delta eta>2.5  "<<cutEff[10]<<"\t:\t"<<((float)cutEff[10]*100.)/(float)cutEff[9]<<std::endl
	  <<"(11) =2 VBF jets   : "<<cutEff[11]<<"\t:\t"<<((float)cutEff[11]*100.)/(float)cutEff[10]<<std::endl
	  <<"(12) Trigger & Events passed all cuts:                       "<<cutEff[12]<<"\t:\t"<<((float)cutEff[12]*100.)/(float)cutEff[11]<<std::endl;//std::cout << "Yield =  " << cutEff[9]*0.00128*WZTree->totalEventWeight<<endl;
  //--------close everything-------------
  delete info; delete gen;
  delete genPartArr; delete electronArr; delete vertexArr;
  delete jetArr;



  // delete vjetArrPuppi; delete vjetAddArrPuppi; 
  delete lheWgtArr;
  outROOT->Write();
  outROOT->Close();
  int t1 = time(NULL);
  printf ("\n==> time to run this code = %0.3f min\n", (float)(t1 - t0)/60.0);
  return(0);
}
