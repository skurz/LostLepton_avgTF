#ifndef Prediction_h
#define Prediction_h

#include "LLTools.h"
#include "SearchBins.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TProfile.h>
#include "TVector2.h"
#include <cmath>
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TKey.h"
#include "TLorentzVector.h"
#include <fstream>
#include "isr/ISRCorrector.h"
#include "btag/BTagCorrector.h"


////////////////////////
//////// Options
//////// SET DEPENDING ON SAMPLE: data, MC, signal!
///////////////////////

// useDeltaPhiCut = 0: no deltaPhiCut
// useDeltaPhiCut = 1: deltaPhiCut
// useDeltaPhiCut = -1: inverted deltaPhiCut
const int useDeltaPhiCut = 1;  //<-check------------------------

const bool runOnData = false;  //<-check------------------------
const bool runOnStandardModelMC = true;  //<-check------------------------
const bool runOnSignalMC = false;  //<-check------------------------

// Use TFs with/without SFs
const bool applySFs = true;

// Use TFs with/without SFs
const double scaleFactorWeight = 35862.351;

// Only needed if running on full nTuples not on Skims (bTag reweighting)
// Does not matter for Data
const bool runOnNtuples = false;
const string path_toSkims("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLe/tree_");

// Useful for T2tt corridor studies
const bool useGenHTMHT = false;

// Do top-pt reweightung
const bool topPTreweight = false;

// pu
const TString path_puHist("pu/PileupHistograms_0721_63mb_pm5.root");
// bTag corrections
const string path_bTagCalib("btag/CSVv2_Moriond17_B_H_mod.csv");
const string path_bTagCalibFastSim("btag/fastsim_csvv2_ttbar_26_1_2017.csv");
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// ISR corrections: NOT RECOMMENDED FOR JAMBOREE -> Might change for Moriond! Just uncomment in Prediction::Init(Tree*) of this file
const TString path_ISRcorr("isr/ISRWeights.root");
// Signal x-sec: "dict_xsec.txt" for gluino pair prod; "dict_xsec_T2.txt" for (anti)stop pair prod.
const string path_xsecT1T5("xsec/dict_xsec.txt");
const string path_xsecT2("xsec/dict_xsec_T2.txt");


////////////////////////
//////// Usually don't have to be changed
///////////////////////

const bool applyFilters=true;


////////////////////////
//////// Don't change anything below
///////////////////////

const double minHT_=300;
//const double minHT_=250;
const double minMHT_=250;
const double minNJets_=1.5;
const double deltaPhi1_=0.5;
const double deltaPhi2_=0.5;
const double deltaPhi3_=0.3;
const double deltaPhi4_=0.3;
const double minDeltaPhiN_=4.0;


class Prediction : public TSelector {
 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  void resetValues();
  bool FiltersPass();

  
  //Options
  bool useTrigger = false;
  bool useTriggerEffWeight = false;
  bool doPUreweighting = false;
  bool doBTagCorr = false;
  bool doISRcorr = false; 
  bool useFilterData = false;

  // output variables
  std::string fname; // for fetching file name
  TString fileName;
  Double_t HTgen_cut = 0;
  TFile* pufile = 0;
  TH1* puhist = 0;

  std::vector<std::pair<double, double>> xsecsT1T5;
  std::vector<std::pair<double, double>> xsecsT2;
  std::vector<std::pair<double, double>> *xsecs = 0;

  // Output
  TH1D* h_Prediction = 0;

  // TFs
  TH1D* h_0L1L_SB = 0;
  TH1D* h_0L1L_SF_SB = 0;

  //open skim file as skimfile
  TH1* h_njetsisr = 0;
  double nEvtsTotal;
  Double_t xsec;
  ISRCorrector *isrcorr = 0;
  TFile* isrfile = 0;
  TH1* h_isr = 0;
  Double_t w_isr;
  Double_t w_pu;
  BTagCorrector *btagcorr = 0;
  std::vector<double> bTagProb;
  std::vector<unsigned int> bTagBins;
  std::vector<unsigned int> bTagBinsQCD;
  Double_t		  topPtSF;
  std::vector<double> topPt;

  TString treeName = " ";


  SearchBins *SearchBins_ =0;
  SearchBins *SearchBinsQCD_ =0;
  SearchBins *SearchBins_BTags_ =0;
  SearchBins *SearchBinsQCD_BTags_ =0;

  Int_t           isoTracksNum;
  UShort_t JetsNum_;
  UShort_t MuonsNoIsoNum_, MuonsNum_;
  UShort_t ElectronsNoIsoNum_, ElectronsNum_;
  Float_t MuonsPt_, MuonsEta_;
  Float_t ElectronsPt_, ElectronsEta_;
  Float_t mtw;
  
  UShort_t Bin_, BinQCD_;

  
  // Declaration of leaf types
  UInt_t          RunNum;
  UInt_t          LumiBlockNum;
  ULong64_t       EvtNum;
  Bool_t           BadChargedCandidateFilter;
  Bool_t           BadPFMuonFilter;
  Int_t           BTags;
  Int_t          CSCTightHaloFilter;
  Double_t        DeltaPhi1;
  Double_t        DeltaPhi2;
  Double_t        DeltaPhi3;
  Double_t        DeltaPhi4;
  Int_t           globalTightHalo2016Filter;
  Int_t           EcalDeadCellTriggerPrimitiveFilter;
  Int_t           eeBadScFilter;
  //Bool_t           eeBadSc4Filter;
  std::vector<TLorentzVector> *GenElectrons=0;
  std::vector<TLorentzVector> *GenMuons=0;
  Int_t          HBHENoiseFilter;
  Int_t          HBHEIsoNoiseFilter;
  Double_t        HT;
  Double_t        GenHT;
  Double_t        GenMHT;
  std::vector<TLorentzVector> *GenJets=0;
  Int_t           isoElectronTracksNum;
  Int_t           isoMuonTracksNum;
  Int_t           isoPionTracksNum;
  Bool_t          JetID;
  std::vector<TLorentzVector> *Jets=0;
  std::vector<double>     *Jets_muonEnergyFraction=0;
  std::vector<double>     *Jets_bDiscriminatorCSV=0;
  std::vector<int>     *Jets_hadronFlavor=0;
  std::vector<int>     *Jets_chargedHadronEnergyFraction=0;
  std::vector<bool>    *Jets_HTMask=0;
  Double_t        METPhi;
  Double_t        MET;
  Double_t        PFCaloMETRatio;
  Double_t        MHT;
  Double_t        MHTPhi;
  Int_t           NJets;
  Int_t           NVtx;
  std::vector<TLorentzVector> *ElectronsNoIso=0;
  std::vector<TLorentzVector> *Electrons=0;
  std::vector<TLorentzVector> *Muons=0;
  std::vector<TLorentzVector> *MuonsNoIso=0;
  std::vector<string>  *TriggerNames=0;
  std::vector<int>    *TriggerPass=0;
  std::vector<int>     *TriggerPrescales=0;
  Double_t        Weight;
  Double_t        puWeight;
  Double_t        madHT;
  Double_t        SusyLSPMass;
  Double_t        SusyMotherMass;
  Double_t        TrueNumInteractions;
  std::vector<double>  *ElectronsNoIso_MT2Activity=0;
  std::vector<double>  *Electrons_MT2Activity=0;
  std::vector<double>  *MuonsNoIso_MT2Activity=0;
  std::vector<double>  *Muons_MT2Activity=0;
  Int_t           NJetsISR;
  vector<TLorentzVector> *GenParticles = 0;
  vector<int>     *GenParticles_PdgId = 0;
  std::vector<double>  *Muons_MTW=0;
  std::vector<double>  *Electrons_MTW=0;
  vector<bool>    *Muons_tightID=0;
  vector<bool>    *Electrons_mediumID=0;
  vector<bool>    *Electrons_tightID=0;



  // List of branches
  TBranch        *b_RunNum=0;   //!
  TBranch        *b_LumiBlockNum=0;   //!
  TBranch        *b_EvtNum=0;   //!
  TBranch        *b_BTags=0;   //!
  TBranch        *b_BadChargedCandidateFilter=0;   //!
  TBranch        *b_BadPFMuonFilter=0;   //!
  TBranch        *b_CSCTightHaloFilter=0;   //!
  TBranch        *b_DeltaPhi1=0;   //!
  TBranch        *b_DeltaPhi2=0;   //!
  TBranch        *b_DeltaPhi3=0;   //!
  TBranch        *b_DeltaPhi4=0;   //!
  TBranch        *b_globalTightHalo2016Filter=0;   //!
  TBranch        *b_EcalDeadCellTriggerPrimitiveFilter=0;   //!
  TBranch        *b_eeBadScFilter=0;   //!
  //TBranch        *b_eeBadSc4Filter=0;   //!
  TBranch        *b_GenElectrons=0;   //!
  TBranch        *b_GenMuons=0;   //!
  TBranch        *b_HBHENoiseFilter=0;   //!
  TBranch        *b_HBHEIsoNoiseFilter=0;   //!
  TBranch        *b_HT=0;   //!
  TBranch        *b_GenHT=0;   //!
  TBranch        *b_GenMHT=0;   //!
  TBranch        *b_GenJets=0;   //!
  TBranch        *b_isoElectronTracksNum=0;   //!
  TBranch        *b_isoMuonTracksNum=0;   //!
  TBranch        *b_isoPionTracksNum=0;   //!
  TBranch        *b_JetID=0;   //!
  TBranch        *b_Jets=0;   //!
  TBranch        *b_Jets_muonEnergyFraction=0;   //!
  TBranch        *b_Jets_bDiscriminatorCSV=0;   //!
  TBranch        *b_Jets_hadronFlavor=0;   //!
  TBranch        *b_Jets_chargedHadronEnergyFraction=0;   //!
  TBranch        *b_Jets_HTMask=0;   //!
  TBranch        *b_METPhi=0;   //!
  TBranch        *b_MET=0;   //!
  TBranch        *b_PFCaloMETRatio=0;   //!
  TBranch        *b_MHT=0;   //!
  TBranch        *b_MHTPhi=0;   //!
  TBranch        *b_NJets=0;   //!
  TBranch        *b_NVtx=0;   //!
  TBranch        *b_ElectronsNoIso=0;   //!
  TBranch        *b_Electrons=0;   //!
  TBranch        *b_Muons=0;   //!
  TBranch        *b_MuonsNoIso=0;   //!
  TBranch        *b_TriggerNames=0;   //!
  TBranch        *b_TriggerPass=0;   //!
  TBranch        *b_TriggerPrescales=0;   //!
  TBranch        *b_Weight=0;   //!
  TBranch        *b_puWeight=0;   //!
  TBranch        *b_madHT=0;
  TBranch        *b_SusyLSPMass=0;
  TBranch        *b_SusyMotherMass=0;
  TBranch        *b_TrueNumInteractions=0;
  TBranch        *b_ElectronsNoIso_MT2Activity=0;   //!
  TBranch        *b_Electrons_MT2Activity=0;   //!
  TBranch        *b_Muons_MT2Activity=0;   //!
  TBranch        *b_MuonsNoIso_MT2Activity=0;   //!
  TBranch        *b_NJetsISR=0;
  TBranch        *b_GenParticles=0;
  TBranch        *b_GenParticles_PdgId=0;
  TBranch        *b_Muons_MTW=0;
  TBranch        *b_Electrons_MTW=0;
  TBranch        *b_Muons_tightID=0;   //!
  TBranch        *b_Electrons_mediumID=0;   //!
  TBranch        *b_Electrons_tightID=0;   //!


  
 Prediction(TTree * /*tree*/ =0) : fChain(0) { }
  virtual ~Prediction() { }
  virtual Int_t   Version() const { return 2; }
  virtual void    Begin(TTree *tree);
  virtual void    SlaveBegin(TTree *tree);
  virtual void    Init(TTree *tree);
  virtual Bool_t  Notify();
  virtual Bool_t  Process(Long64_t entry);
  virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
  virtual void    SetOption(const char *option) { fOption = option; }
  virtual void    SetObject(TObject *obj) { fObject = obj; }
  virtual void    SetInputList(TList *input) { fInput = input; }
  virtual TList  *GetOutputList() const { return fOutput; }
  virtual void    SlaveTerminate();
  virtual void    Terminate();
  
  ClassDef(Prediction,0);
};

#endif

#ifdef Prediction_cxx
void Prediction::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).
  
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fChain->SetMakeClass(1);


  ////////////////////////
  //////// Options
  ///////////////////////

  if(!((runOnData && !runOnStandardModelMC && !runOnSignalMC) || (!runOnData && runOnStandardModelMC && !runOnSignalMC) || (!runOnData && !runOnStandardModelMC && runOnSignalMC))){
    fprintf(stderr, "CHECK OPTIONS! EITHER RUN ON DATA, MC, OR SIGNAL!!!");
    return;
  }

  // Seta data specific options, e.g. Trigger
  // Apply trigger
  if(runOnData) useTrigger = true;
  // Apply weights if trigger not simulated
  //if(runOnStandardModelMC) useTriggerEffWeight = true; // not derived yet
  if(runOnSignalMC && !useGenHTMHT) useTriggerEffWeight = true;
  // Do PU reweighting. true for signal scan
  if(runOnSignalMC) doPUreweighting = true;
  //if(runOnStandardModelMC) doPUreweighting = true;
  // bTag corrections. Use for signal scan
  if(!runOnData) doBTagCorr = true;
  // ISR corrections.
  if(runOnSignalMC) doISRcorr = true; //<-check---------------------------------------

  // useFilterData = true; unless you want to run without MET filters
  // useFilterData = false; For FastSim Samples, e.g. Signal Scans! Met filters not simulated
  if(runOnStandardModelMC || runOnData) useFilterData = true;



  ////////////////////////
  //////// End Options
  ///////////////////////

  // Open histograms with TFs
  TFile *TF_histFile = TFile::Open("TF.root", "READ");
  h_0L1L_SB = (TH1D*) TF_histFile->Get("h_0L1L_SB")->Clone();
  h_0L1L_SF_SB = (TH1D*) TF_histFile->Get("h_0L1L_SF_SB")->Clone();

  TChain* temp = (TChain*)fChain;
  std::string infname=temp->GetFile()->GetName();

  TFile* skimfile = temp->GetFile();

  std::string baseName(infname);
  size_t pos=baseName.rfind("/");
  if(pos!=std::string::npos){
    if(pos!=baseName.size()-1){
      baseName.erase(0,pos+1);
    }
  }
  pos=baseName.rfind(".root");
  if(pos!=std::string::npos){
    if(pos!=baseName.size()-1){
      baseName.erase(pos);
    }
  }
 
  fname = baseName+"_Exp.root";

  TString option = GetOption();
  TObjArray *optionArray = option.Tokenize(",");
  fileName = fname.c_str();

  TString fileNameString = "";
  TString HTcutString = "";

  if(!optionArray->IsEmpty()){
    fileNameString = ((TObjString *)(optionArray->At(0)))->String();
    if(optionArray->GetEntries() > 1) HTcutString = ((TObjString *)(optionArray->At(1)))->String();
  }

  fileNameString = fileNameString.Strip(TString::kBoth, ' ').String();
  HTcutString = HTcutString.Strip(TString::kBoth, ' ').String();

  if(fileNameString!="*" && fileNameString!="") fileName = fileNameString;
  if(HTcutString!="" && HTcutString.IsFloat()) HTgen_cut = HTcutString.Atof();

  std::cout << "madHT cut: " << HTgen_cut << std::endl;

  std::cout << "Saving file to: " << fileName << std::endl;

  std::ifstream signal_xsecT1T5(path_xsecT1T5);
  std::string str_xsecT1T5;
  while (std::getline(signal_xsecT1T5, str_xsecT1T5))
  {
    TObjArray *tokens = TString(str_xsecT1T5).Tokenize(",");
    //std::cout<<((TObjString *)(tokens->At(0)))->String()<<"; "<<((TObjString *)(tokens->At(1)))->String()<<";"<<std::endl;
    xsecsT1T5.push_back(std::make_pair(std::atof(((TObjString *)(tokens->At(0)))->String()), std::atof(((TObjString *)(tokens->At(1)))->String())));
  }

  std::ifstream signal_xsecT2(path_xsecT2);
  std::string str_xsecT2;
  while (std::getline(signal_xsecT2, str_xsecT2))
  {
    TObjArray *tokens = TString(str_xsecT2).Tokenize(",");
    //std::cout<<((TObjString *)(tokens->At(0)))->String()<<"; "<<((TObjString *)(tokens->At(1)))->String()<<";"<<std::endl;
    xsecsT2.push_back(std::make_pair(std::atof(((TObjString *)(tokens->At(0)))->String()), std::atof(((TObjString *)(tokens->At(1)))->String())));
  }

  if(runOnSignalMC){
    // ISR setup
    isrfile = TFile::Open(path_ISRcorr, "READ");
    h_isr = (TH1*)isrfile->Get("isr_weights_central");
  }

  if(doPUreweighting){
    pufile = TFile::Open(path_puHist,"READ");
    puhist = (TH1*)pufile->Get("pu_weights_central");
  }

  fChain->SetBranchStatus("*",0);

  fChain->SetBranchStatus("RunNum", 1);
  fChain->SetBranchStatus("LumiBlockNum", 1);
  fChain->SetBranchStatus("EvtNum", 1);
  fChain->SetBranchStatus("BTags", 1);
  fChain->SetBranchStatus("DeltaPhi1", 1);
  fChain->SetBranchStatus("DeltaPhi2", 1);
  fChain->SetBranchStatus("DeltaPhi3", 1);
  fChain->SetBranchStatus("DeltaPhi4", 1);
  if(!runOnSignalMC){
    //fChain->SetBranchStatus("CSCTightHaloFilter", 1);
    fChain->SetBranchStatus("EcalDeadCellTriggerPrimitiveFilter", 1);
    fChain->SetBranchStatus("eeBadScFilter", 1);
    //fChain->SetBranchStatus("eeBadSc4Filter", 1);    
    fChain->SetBranchStatus("HBHENoiseFilter", 1);
    fChain->SetBranchStatus("HBHEIsoNoiseFilter", 1);
    if(runOnData){
      fChain->SetBranchStatus("globalTightHalo2016Filter", 1);
      fChain->SetBranchStatus("BadChargedCandidateFilter", 1);
      fChain->SetBranchStatus("BadPFMuonFilter", 1);
    }
  }
  fChain->SetBranchStatus("Electrons", 1);
  fChain->SetBranchStatus("HT", 1);
  fChain->SetBranchStatus("isoElectronTracks", 1);
  fChain->SetBranchStatus("isoMuonTracks", 1);
  fChain->SetBranchStatus("isoPionTracks", 1);
  fChain->SetBranchStatus("JetID", 1);
  fChain->SetBranchStatus("Jets", 1);
  fChain->SetBranchStatus("Jets_HTMask", 1);
  fChain->SetBranchStatus("METPhi", 1);
  fChain->SetBranchStatus("MET", 1);
  fChain->SetBranchStatus("PFCaloMETRatio", 1);
  fChain->SetBranchStatus("MHT", 1);
  fChain->SetBranchStatus("MHTPhi", 1);
  fChain->SetBranchStatus("Muons", 1);
  fChain->SetBranchStatus("NJets", 1);
  fChain->SetBranchStatus("NVtx", 1);
  fChain->SetBranchStatus("ElectronsNoIso", 1);
  fChain->SetBranchStatus("Electrons", 1);
  fChain->SetBranchStatus("Muons", 1);
  fChain->SetBranchStatus("MuonsNoIso", 1);
  fChain->SetBranchStatus("TriggerNames", 1);
  fChain->SetBranchStatus("TriggerPass", 1);
  fChain->SetBranchStatus("TriggerPrescales", 1);
  fChain->SetBranchStatus("Jets_muonEnergyFraction", 1);
  fChain->SetBranchStatus("Jets_bDiscriminatorCSV", 1);
  if(topPTreweight){
    fChain->SetBranchStatus("GenParticles", 1);
    fChain->SetBranchStatus("GenParticles_PdgId", 1);
    fChain->SetBranchStatus("GenElectrons", 1);
    fChain->SetBranchStatus("GenMuons", 1);
  }  
  

  if(!runOnData){
    fChain->SetBranchStatus("Weight", 1);
    fChain->SetBranchStatus("Jets_hadronFlavor", 1);
    fChain->SetBranchStatus("madHT", 1);
    fChain->SetBranchStatus("TrueNumInteractions", 1);
  } 
  if(runOnSignalMC){
    fChain->SetBranchStatus("SusyLSPMass", 1);
    fChain->SetBranchStatus("SusyMotherMass", 1);
    fChain->SetBranchStatus("NJetsISR", 1);
    fChain->SetBranchStatus("GenJets", 1);
    fChain->SetBranchStatus("Jets_chargedHadronEnergyFraction", 1);
  }

  //if(useGenHTMHT){
  //  fChain->SetBranchStatus("GenHT", 1);
  //  fChain->SetBranchStatus("GenMHT", 1);
  //}

  fChain->SetBranchStatus("ElectronsNoIso_MT2Activity",1);
  fChain->SetBranchStatus("Electrons_MT2Activity", 1);
  fChain->SetBranchStatus("Muons_MT2Activity",1);
  fChain->SetBranchStatus("MuonsNoIso_MT2Activity", 1);
  fChain->SetBranchStatus("Muons_MTW", 1);
  fChain->SetBranchStatus("Electrons_MTW", 1);
  fChain->SetBranchStatus("Muons_tightID", 1);
  fChain->SetBranchStatus("Electrons_mediumID", 1);
  fChain->SetBranchStatus("Electrons_tightID", 1);

  fChain->SetBranchAddress("RunNum", &RunNum, &b_RunNum);
  fChain->SetBranchAddress("LumiBlockNum", &LumiBlockNum, &b_LumiBlockNum);
  fChain->SetBranchAddress("EvtNum", &EvtNum, &b_EvtNum);
  fChain->SetBranchAddress("BTags", &BTags, &b_BTags);
  fChain->SetBranchAddress("DeltaPhi1", &DeltaPhi1, &b_DeltaPhi1);
  fChain->SetBranchAddress("DeltaPhi2", &DeltaPhi2, &b_DeltaPhi2);
  fChain->SetBranchAddress("DeltaPhi3", &DeltaPhi3, &b_DeltaPhi3);
  fChain->SetBranchAddress("DeltaPhi4", &DeltaPhi4, &b_DeltaPhi4);
  if(!runOnSignalMC){
    //fChain->SetBranchAddress("CSCTightHaloFilter", &CSCTightHaloFilter, &b_CSCTightHaloFilter);
    fChain->SetBranchAddress("EcalDeadCellTriggerPrimitiveFilter", &EcalDeadCellTriggerPrimitiveFilter, &b_EcalDeadCellTriggerPrimitiveFilter);
    fChain->SetBranchAddress("eeBadScFilter", &eeBadScFilter, &b_eeBadScFilter);
    //fChain->SetBranchAddress("eeBadSc4Filter", &eeBadSc4Filter, &b_eeBadSc4Filter);    
    fChain->SetBranchAddress("HBHENoiseFilter", &HBHENoiseFilter, &b_HBHENoiseFilter);
    fChain->SetBranchAddress("HBHEIsoNoiseFilter", &HBHEIsoNoiseFilter, &b_HBHEIsoNoiseFilter);
    if(runOnData){
      fChain->SetBranchAddress("globalTightHalo2016Filter", &globalTightHalo2016Filter, &b_globalTightHalo2016Filter);
      fChain->SetBranchAddress("BadChargedCandidateFilter", &BadChargedCandidateFilter, &b_BadChargedCandidateFilter);
      fChain->SetBranchAddress("BadPFMuonFilter", &BadPFMuonFilter, &b_BadPFMuonFilter);
    }
  }
  fChain->SetBranchAddress("HT", &HT, &b_HT);
  fChain->SetBranchAddress("isoElectronTracks", &isoElectronTracksNum, &b_isoElectronTracksNum);
  fChain->SetBranchAddress("isoMuonTracks", &isoMuonTracksNum, &b_isoMuonTracksNum);
  fChain->SetBranchAddress("isoPionTracks", &isoPionTracksNum, &b_isoPionTracksNum);
  fChain->SetBranchAddress("JetID", &JetID, &b_JetID);
  fChain->SetBranchAddress("Jets", &Jets, &b_Jets);
  fChain->SetBranchAddress("Jets_HTMask", &Jets_HTMask, &b_Jets_HTMask);
  fChain->SetBranchAddress("METPhi", &METPhi, &b_METPhi);
  fChain->SetBranchAddress("MET", &MET, &b_MET);
  fChain->SetBranchAddress("PFCaloMETRatio", &PFCaloMETRatio, &b_PFCaloMETRatio);
  fChain->SetBranchAddress("MHT", &MHT, &b_MHT);
  fChain->SetBranchAddress("MHTPhi", &MHTPhi, &b_MHTPhi);
  fChain->SetBranchAddress("NJets", &NJets, &b_NJets);
  fChain->SetBranchAddress("NVtx", &NVtx, &b_NVtx);
  fChain->SetBranchAddress("ElectronsNoIso", &ElectronsNoIso, &b_ElectronsNoIso);
  fChain->SetBranchAddress("Electrons", &Electrons, &b_Electrons);
  fChain->SetBranchAddress("Muons", &Muons, &b_Muons);
  fChain->SetBranchAddress("MuonsNoIso", &MuonsNoIso, &b_MuonsNoIso);
  fChain->SetBranchAddress("TriggerNames", &TriggerNames, &b_TriggerNames);
  fChain->SetBranchAddress("TriggerPass", &TriggerPass, &b_TriggerPass);
  fChain->SetBranchAddress("TriggerPrescales", &TriggerPrescales, &b_TriggerPrescales);
  fChain->SetBranchAddress("Jets_muonEnergyFraction", &Jets_muonEnergyFraction, &b_Jets_muonEnergyFraction);
  fChain->SetBranchAddress("Jets_bDiscriminatorCSV", &Jets_bDiscriminatorCSV, &b_Jets_bDiscriminatorCSV);
  if(topPTreweight){
    fChain->SetBranchAddress("GenParticles", &GenParticles, &b_GenParticles);
    fChain->SetBranchAddress("GenParticles_PdgId", &GenParticles_PdgId, &b_GenParticles_PdgId);
    fChain->SetBranchAddress("GenElectrons", &GenElectrons, &b_GenElectrons);
    fChain->SetBranchAddress("GenMuons", &GenMuons, &b_GenMuons);
  }  
  if(!runOnData){
    fChain->SetBranchAddress("Weight", &Weight, &b_Weight);
    fChain->SetBranchAddress("Jets_hadronFlavor", &Jets_hadronFlavor, &b_Jets_hadronFlavor);
    fChain->SetBranchAddress("madHT", &madHT, &b_madHT);
    fChain->SetBranchAddress("TrueNumInteractions", &TrueNumInteractions, &b_TrueNumInteractions);
  }
  if(runOnSignalMC){
    fChain->SetBranchAddress("SusyLSPMass", &SusyLSPMass, &b_SusyLSPMass);
    fChain->SetBranchAddress("SusyMotherMass", &SusyMotherMass, &b_SusyMotherMass);
    fChain->SetBranchAddress("NJetsISR", &NJetsISR, &b_NJetsISR);
    fChain->SetBranchAddress("GenJets", &GenJets, &b_GenJets);
    fChain->SetBranchAddress("Jets_chargedHadronEnergyFraction", &Jets_chargedHadronEnergyFraction, &b_Jets_chargedHadronEnergyFraction);
  }

  //if(useGenHTMHT){
  //  fChain->SetBranchAddress("GenHT", &GenHT, &b_GenHT);
  //  fChain->SetBranchAddress("GenMHT", &GenMHT, &b_GenMHT);
  //}

  fChain->SetBranchAddress("ElectronsNoIso_MT2Activity", &ElectronsNoIso_MT2Activity, &b_ElectronsNoIso_MT2Activity);
  fChain->SetBranchAddress("Electrons_MT2Activity", &Electrons_MT2Activity, &b_Electrons_MT2Activity);
  fChain->SetBranchAddress("Muons_MT2Activity", &Muons_MT2Activity, &b_Muons_MT2Activity);
  fChain->SetBranchAddress("MuonsNoIso_MT2Activity", &MuonsNoIso_MT2Activity, &b_MuonsNoIso_MT2Activity);
  fChain->SetBranchAddress("Muons_MTW", &Muons_MTW, &b_Muons_MTW);
  fChain->SetBranchAddress("Electrons_MTW", &Electrons_MTW, &b_Electrons_MTW);
  fChain->SetBranchAddress("Muons_tightID", &Muons_tightID, &b_Muons_tightID);
  fChain->SetBranchAddress("Electrons_mediumID", &Electrons_mediumID, &b_Electrons_mediumID);
  fChain->SetBranchAddress("Electrons_tightID", &Electrons_tightID, &b_Electrons_tightID);
}


Bool_t Prediction::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.
  
  return kTRUE;
}

#endif // #ifdef Prediction_cxx
