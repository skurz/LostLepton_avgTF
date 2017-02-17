#include <TChain.h>
#include "TProofServ.h"
#include "TProof.h"
#include <iostream>


using std::vector;

//needed to write vector<TLorentzVector> to tree
#ifdef __CINT__
#pragma link C++ class std::vector<TLorentzVector>+;
#endif

void MakePrediction_Data()
{
	gSystem->Load("libPhysics.so");
  	gInterpreter->GenerateDictionary("vector<TLorentzVector>","TLorentzVector.h;vector");

	bool useTProof = false;

	TChain *Effchain = new TChain("tree");
	Effchain->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLm/tree_MET_2016B.root");	
	Effchain->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLe/tree_MET_2016B.root");
	Effchain->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLm/tree_MET_2016C.root");	
	Effchain->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLe/tree_MET_2016C.root");
	Effchain->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLm/tree_MET_2016D.root");	
	Effchain->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLe/tree_MET_2016D.root");
	Effchain->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLm/tree_MET_2016E.root");	
	Effchain->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLe/tree_MET_2016E.root");
	Effchain->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLm/tree_MET_2016F.root");	
	Effchain->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLe/tree_MET_2016F.root");
	Effchain->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLm/tree_MET_2016G.root");	
	Effchain->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLe/tree_MET_2016G.root");
	Effchain->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLm/tree_MET_2016H2.root");	
	Effchain->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLe/tree_MET_2016H2.root");
	Effchain->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLm/tree_MET_2016H3.root");	
	Effchain->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLe/tree_MET_2016H3.root");

//	TChain *Effchain = new TChain("TreeMaker2/PreSelection");
//	Effchain->Add("/nfs/dust/cms/user/kurzsimo/LostLepton/data_v12/Run2016*.root");

	TProof *proof = NULL;

    if(useTProof){
    	proof = TProof::Open("workers=10");
    	Effchain->SetProof();     	
    }

	Effchain->Process("Prediction.C++g", "Prediction_data.root");

    delete proof;
}
