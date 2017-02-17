#include <TChain.h>
#include <vector>
#include "TProofServ.h"
#include "TProof.h"
#include "TLorentzVector.h"
#include "TROOT.h"
#include <iostream>


using std::vector;

//needed to write vector<TLorentzVector> to tree
#ifdef __CINT__
#pragma link C++ class std::vector<TLorentzVector>+;
#endif


void MakePrediction_separate()
{
  gSystem->Load("libPhysics.so");
  gInterpreter->GenerateDictionary("vector<TLorentzVector>","TLorentzVector.h;vector");

  gROOT->ProcessLine(".L Prediction.C+");
  
  const int nChains = 33;
  TChain *Effchain[nChains];
  for(Int_t i=0; i<nChains; i++){
    Effchain[i] = new TChain("tree");
  }

  Effchain[0]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLe/tree_TTJets_SingleLeptFromT.root");
  Effchain[0]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLm/tree_TTJets_SingleLeptFromT.root");
  Effchain[1]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLe/tree_TTJets_SingleLeptFromTbar.root");
  Effchain[1]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLm/tree_TTJets_SingleLeptFromTbar.root");
  Effchain[2]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLe/tree_TTJets_DiLept.root");
  Effchain[2]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLm/tree_TTJets_DiLept.root");

  Effchain[3]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLe/tree_TTJets_HT-600to800.root");
  Effchain[3]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLm/tree_TTJets_HT-600to800.root");
  Effchain[4]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLe/tree_TTJets_HT-800to1200.root");
  Effchain[4]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLm/tree_TTJets_HT-800to1200.root");
  Effchain[5]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLe/tree_TTJets_HT-1200to2500.root");
  Effchain[5]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLm/tree_TTJets_HT-1200to2500.root");
  Effchain[6]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLe/tree_TTJets_HT-2500toInf.root");
  Effchain[6]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLm/tree_TTJets_HT-2500toInf.root");

  Effchain[7]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLe/tree_WJetsToLNu_HT-100to200.root");
  Effchain[7]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLm/tree_WJetsToLNu_HT-100to200.root");
  Effchain[8]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLe/tree_WJetsToLNu_HT-200to400.root");
  Effchain[8]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLm/tree_WJetsToLNu_HT-200to400.root");
  Effchain[9]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLe/tree_WJetsToLNu_HT-400to600.root");
  Effchain[9]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLm/tree_WJetsToLNu_HT-400to600.root");
  Effchain[10]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLe/tree_WJetsToLNu_HT-600to800.root");
  Effchain[10]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLm/tree_WJetsToLNu_HT-600to800.root");
  Effchain[11]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLe/tree_WJetsToLNu_HT-800to1200.root");
  Effchain[11]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLm/tree_WJetsToLNu_HT-800to1200.root");
  Effchain[12]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLe/tree_WJetsToLNu_HT-1200to2500.root");
  Effchain[12]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLm/tree_WJetsToLNu_HT-1200to2500.root");
  Effchain[13]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLe/tree_WJetsToLNu_HT-2500toInf.root");
  Effchain[13]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLm/tree_WJetsToLNu_HT-2500toInf.root");

  Effchain[14]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLe/tree_ST_s-channel.root");
  Effchain[14]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLm/tree_ST_s-channel.root");
  Effchain[15]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLe/tree_ST_t-channel_antitop.root");
  Effchain[15]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLm/tree_ST_t-channel_antitop.root");
  Effchain[16]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLe/tree_ST_t-channel_top.root");
  Effchain[16]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLm/tree_ST_t-channel_top.root");
  Effchain[17]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLe/tree_ST_tW_antitop.root");
  Effchain[17]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLm/tree_ST_tW_antitop.root");
  Effchain[18]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLe/tree_ST_tW_top.root");
  Effchain[18]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLm/tree_ST_tW_top.root");

  Effchain[19]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLe/tree_TTGJets.root");
  Effchain[19]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLm/tree_TTGJets.root");
  Effchain[20]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLe/tree_TTTT.root");
  Effchain[20]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLm/tree_TTTT.root");
  Effchain[21]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLe/tree_TTWJetsToLNu.root");
  Effchain[21]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLm/tree_TTWJetsToLNu.root");
  Effchain[22]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLe/tree_TTWJetsToQQ.root");
  Effchain[22]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLm/tree_TTWJetsToQQ.root");
  Effchain[23]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLe/tree_TTZToLLNuNu.root");
  Effchain[23]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLm/tree_TTZToLLNuNu.root");
  Effchain[24]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLe/tree_TTZToQQ.root");
  Effchain[24]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLm/tree_TTZToQQ.root");
  Effchain[25]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLe/tree_WWTo1L1Nu2Q.root");
  Effchain[25]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLm/tree_WWTo1L1Nu2Q.root");
  Effchain[26]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLe/tree_WWTo2L2Nu.root");
  Effchain[26]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLm/tree_WWTo2L2Nu.root");
  Effchain[27]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLe/tree_WWZ.root");
  Effchain[27]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLm/tree_WWZ.root");
  Effchain[28]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLe/tree_WZTo1L1Nu2Q.root");
  Effchain[28]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLm/tree_WZTo1L1Nu2Q.root");
  Effchain[29]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLe/tree_WZTo1L3Nu.root");
  Effchain[29]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLm/tree_WZTo1L3Nu.root");
  Effchain[30]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLe/tree_WZZ.root");
  Effchain[30]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLm/tree_WZZ.root");
  Effchain[31]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLe/tree_ZZTo2L2Q.root");
  Effchain[31]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLm/tree_ZZTo2L2Q.root");
  Effchain[32]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLe/tree_ZZZ.root");
  Effchain[32]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLm/tree_ZZZ.root");
  //Effchain[33]->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLe/tree_ZZTo2Q2Nu.root");


  for(Int_t i=0; i<nChains; i++){ //i<nChains i>2
    std::cout<<"Processing Tree: "<<i<<std::endl;
    if(i>2)
        Effchain[i]->Process("Prediction", TString::Format("Prediction_separate/Prediction_%d.root",i));
    else    Effchain[i]->Process("Prediction", TString::Format("Prediction_separate/Prediction_%d.root, 600",i)); // folder already has to exist for the time beeing...
    
  }
}
