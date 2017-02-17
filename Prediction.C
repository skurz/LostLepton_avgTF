#define Prediction_cxx

#include "Prediction.h"
#include <TH2.h>
#include <TStyle.h>
#include <iostream>


void Prediction::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

}

void Prediction::SlaveBegin(TTree * /*tree*/)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).
  

  SearchBins_ = new SearchBins(false);
  SearchBinsQCD_ = new SearchBins(true);
  SearchBins_BTags_ = new SearchBins(false);
  SearchBinsQCD_BTags_ = new SearchBins(true);

  bTagBins = {0, 0, 0, 0};
  bTagBinsQCD = {0, 0, 0, 0};

  // Initialize Histograms
  TH1::SetDefaultSumw2();
  unsigned nSB = SearchBins_->GetNbins();
  h_Prediction = new TH1D("h_Prediction", "h_Prediction", nSB, 0.5, nSB+0.5);

  GetOutputList()->Add(h_Prediction);

  std::cout<<"Run on Data: "<<runOnData<<std::endl;
  std::cout<<"Run on SM MC: "<<runOnStandardModelMC<<std::endl;
  std::cout<<"Run on Signal MC: "<<runOnSignalMC<<std::endl;
  std::cout<<"----------------"<<std::endl;
  std::cout<<"DeltaPhi Cut: "<<useDeltaPhiCut<<std::endl;
  std::cout<<"----------------"<<std::endl;
}

Bool_t Prediction::Process(Long64_t entry)
{
  resetValues();
  fChain->GetTree()->GetEntry(entry);

  if(HTgen_cut > 0.01) if(madHT > HTgen_cut) return kTRUE;

  MuonsNum_ = Muons->size();
  ElectronsNum_ = Electrons->size();

  if((MuonsNum_+ElectronsNum_) !=1) return kTRUE;
  
  if(HT<minHT_ || MHT< minMHT_ || NJets < minNJets_  ) return kTRUE;
  if(useDeltaPhiCut == 1) if(DeltaPhi1 < deltaPhi1_ || DeltaPhi2 < deltaPhi2_ || DeltaPhi3 < deltaPhi3_ || DeltaPhi4 < deltaPhi4_) return kTRUE;
  if(useDeltaPhiCut == -1) if(!(DeltaPhi1 < deltaPhi1_ || DeltaPhi2 < deltaPhi2_ || DeltaPhi3 < deltaPhi3_ || DeltaPhi4 < deltaPhi4_)) return kTRUE;
  if(applyFilters &&  !FiltersPass() ) return kTRUE;

  if(MuonsNum_==1 && ElectronsNum_==0){
    mtw =  Muons_MTW->at(0);
  }else if(MuonsNum_==0 && ElectronsNum_==1){
    mtw =  Electrons_MTW->at(0);
  }
  if(mtw > 100) return kTRUE;

  isoTracksNum = isoElectronTracksNum + isoMuonTracksNum + isoPionTracksNum;

  // Signal triggers
  if(useTrigger) if(!TriggerPass->at(42) && !TriggerPass->at(43) &&!TriggerPass->at(44) && !TriggerPass->at(46) && !TriggerPass->at(47) && !TriggerPass->at(48)) return kTRUE;

  Bin_ = SearchBins_->GetBinNumber(HT,MHT,NJets,BTags);
  BinQCD_ = SearchBinsQCD_->GetBinNumber(HT,MHT,NJets,BTags);

  if(Bin_ > 900 && BinQCD_ > 900) return kTRUE;

  bTagProb = {1, 0, 0, 0};
  bTagBins = {Bin_, 0, 0, 0};
  bTagBinsQCD = {BinQCD_, 0, 0, 0};


  if(topPTreweight){
    if(GenParticles->size() != GenParticles_PdgId->size()){
      std::cout << "Cannot do top-pT reweighting!"<< std::endl; 
    }else{
      for(unsigned iGen = 0; iGen < GenParticles->size(); iGen++){
        if(std::abs(GenParticles_PdgId->at(iGen)) == 6){
          topPt.push_back(GenParticles->at(iGen).Pt());
        }
      }
      
      // https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting#Example
      if(topPt.size() == 2){
        // dilept
        if(GenElectrons->size() + GenMuons->size() == 2){
          topPtSF = std::sqrt(std::exp(0.148-0.00129*(topPt.at(0) < 400. ? topPt.at(0) : 400.))*std::exp(0.148-0.00129*(topPt.at(1) < 400. ? topPt.at(1) : 400.)));
        // singlelept
        }else if(GenElectrons->size() + GenMuons->size() == 1){
          topPtSF = std::sqrt(std::exp(0.159-0.00141*(topPt.at(0) < 400. ? topPt.at(0) : 400.))*std::exp(0.159-0.00141*(topPt.at(1) < 400. ? topPt.at(1) : 400.)));
        //had
        }else{
          // Usually non-promt (in hadTau evts): use average SF
          topPtSF = std::sqrt(std::exp(0.156-0.00137*(topPt.at(0) < 400. ? topPt.at(0) : 400.))*std::exp(0.156-0.00137*(topPt.at(1) < 400. ? topPt.at(1) : 400.)));
          //std::cout << "Cannot do top-pT reweighting! No leptonic top found."<< std::endl; 
        }
      }else{
        topPtSF = -1;
        std::cout << "Cannot do top-pT reweighting! More/Less than 2 tops found."<< std::endl; 
      }

    }

    // Normalization tested on SingleLept and DiLept samples (varies from ~98.9x-99.0x)
        topPtSF /= 0.99;
  }

  if(!runOnData){
    TString currentTree = TString(fChain->GetCurrentFile()->GetName());

    if(currentTree != treeName){
      treeName = currentTree;

      if(doISRcorr){
        h_njetsisr = (TH1*) fChain->GetCurrentFile()->Get("NJetsISR");
        if(isrcorr!=0){
          delete isrcorr;
          isrcorr = 0;
        }
        isrcorr = new ISRCorrector();
        isrcorr->SetWeights(h_isr,h_njetsisr);
      }

      if(doBTagCorr){
        if(btagcorr!=0){
          delete btagcorr;
          btagcorr = 0;
        }
        btagcorr = new BTagCorrector();

        if(!runOnNtuples) btagcorr->SetEffs(fChain->GetCurrentFile());
        else{
          TObjArray *optionArray = currentTree.Tokenize("/");
          TString currFileName = ((TObjString *)(optionArray->At(optionArray->GetEntries()-1)))->String();
          TFile *skimFile = TFile::Open(path_toSkims+currFileName, "READ");
          btagcorr->SetEffs(skimFile);
        }

        btagcorr->SetCalib(path_bTagCalib);        
        if(runOnSignalMC){
          btagcorr->SetCalibFastSim(path_bTagCalibFastSim);
          btagcorr->SetFastSim(true);
        }
        else btagcorr->SetFastSim(false);
      }

      if(runOnSignalMC){
        if((std::string(currentTree.Data()).find(std::string("T1"))) != std::string::npos || (std::string(currentTree.Data()).find(std::string("T5"))) != std::string::npos){
          xsecs = &xsecsT1T5;
          //std::cout<<"Using xsecs for gluino pair production!"<<std::endl;
        }else if((std::string(currentTree.Data()).find(std::string("T2"))) != std::string::npos){
          xsecs = &xsecsT2;
          //std::cout<<"Using xsecs for stop pair production!"<<std::endl;
        }else{
          std::cout<<"No valid dictionary with xsecs found!"<<std::endl;
          return kTRUE;
        }
      }
    }

    if(runOnSignalMC){
      TH1F *nEventProc = (TH1F*)fChain->GetCurrentFile()->Get("nEventProc");
      TH1F *nEventNeg = (TH1F*)fChain->GetCurrentFile()->Get("nEventNeg");
      nEvtsTotal = nEventProc->GetBinContent(1) - 2*nEventNeg->GetBinContent(1);

      xsec = 0;
      for (std::vector<std::pair<double, double>>::iterator it = xsecs->begin() ; it != xsecs->end(); ++it){
        if(std::abs(SusyMotherMass - it->first) < 0.1){
          xsec = it->second;
          break;
        }
      }
      if(xsec < 1e-10){
        std::cout<<"No valid xsec found!"<<std::endl;
        return kTRUE;
      }

      Weight = xsec / nEvtsTotal;
      if(Weight < 0) Weight *= -1;
    }   

    if(doISRcorr){
      w_isr = isrcorr->GetCorrection(NJetsISR);
      Weight *= w_isr;
    }

    if(doBTagCorr){
      bTagProb = btagcorr->GetCorrections(Jets,Jets_hadronFlavor,Jets_HTMask);
      bTagBins = {SearchBins_BTags_->GetBinNumber(HT,MHT,NJets,0), SearchBins_BTags_->GetBinNumber(HT,MHT,NJets,1), SearchBins_BTags_->GetBinNumber(HT,MHT,NJets,2), NJets < 3 ? 999 : SearchBins_BTags_->GetBinNumber(HT,MHT,NJets,3)};  
      bTagBinsQCD = {SearchBinsQCD_BTags_->GetBinNumber(HT,MHT,NJets,0), SearchBinsQCD_BTags_->GetBinNumber(HT,MHT,NJets,1), SearchBinsQCD_BTags_->GetBinNumber(HT,MHT,NJets,2), NJets < 3 ? 999 : SearchBinsQCD_BTags_->GetBinNumber(HT,MHT,NJets,3)};  
    }
  }

  if(runOnSignalMC){
    //Account for efficiency of JetID since we cannot apply it on fastSim
    Weight *= 0.99;
  }

  if(useTriggerEffWeight){
    if(runOnSignalMC){
      Weight *= GetSignalTriggerEffWeight(MHT);
    }else{
      Weight *= GetTriggerEffWeight(MHT);
    }
  }

  if(doPUreweighting){
    w_pu = puhist->GetBinContent(puhist->GetXaxis()->FindBin(min(TrueNumInteractions,puhist->GetBinLowEdge(puhist->GetNbinsX()+1))));
    Weight *= w_pu;
  }

  if(runOnData) Weight = 1.;
  else Weight *= scaleFactorWeight;


  int nLoops = 1;
  if(doBTagCorr) nLoops = (NJets == 2 ? 3 : 4);
  for(int i = 0; i < nLoops; i++){
    double WeightBtagProb = Weight*bTagProb.at(i);
    unsigned bTagBin = bTagBins.at(i);
    unsigned bTagBinQCD = bTagBinsQCD.at(i);

    double TF = -1;
    if(applySFs){
      TF = h_0L1L_SF_SB->GetBinContent(bTagBinQCD);
      if(TF < 0) TF = h_0L1L_SB->GetBinContent(bTagBinQCD);
    }else{
      TF = h_0L1L_SB->GetBinContent(bTagBinQCD);
    }

    h_Prediction->Fill(bTagBin, WeightBtagProb*TF);
  }


  return kTRUE;
}

void Prediction::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

  //std::cout<<"--- QCD binning ---"<<std::endl;
  //SearchBinsQCD_->PrintUsed();

  std::cout<<"--- Search bins ---"<<std::endl;
  SearchBins_->PrintUsed();  
}

void Prediction::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.


  h_Prediction = dynamic_cast<TH1D*>(GetOutputList()->FindObject("h_Prediction"));

  TFile *outPutFile = new TFile(fileName,"RECREATE"); ;
  outPutFile->cd();

  h_Prediction->Write();

  outPutFile->Close();

  cout << "Saved output to " << fileName << endl;
  
}

void Prediction::resetValues()
{
    mtw = 0;
}

bool Prediction::FiltersPass()
{
  bool result=true;
  if(useFilterData){
    if(HBHENoiseFilter!=1) result=false;
    if(HBHEIsoNoiseFilter!=1) result=false;
    if(EcalDeadCellTriggerPrimitiveFilter!=1) result=false;    
    if(eeBadScFilter!=1) result=false;

    if(runOnData){
      if(!BadChargedCandidateFilter) result=false;
      if(!BadPFMuonFilter) result=false;
      if(globalTightHalo2016Filter!=1) result=false;
    }    
  }
  if(NVtx<=0) result=false;

  // Do not apply on fastSim samples!
  if(!runOnSignalMC) if(!JetID) result=false;

  // Preliminary filters
  if(PFCaloMETRatio>5) result=false;

  // Check efficiency of filter
  
  if(result)
  for(unsigned j = 0; j < Jets->size(); j++){
    if(TMath::IsNaN(Jets->at(j).Phi()-METPhi)) result=false;
    if(Jets->at(j).Pt()>200 && Jets_muonEnergyFraction->at(j)>0.5 && (TVector2::Phi_mpi_pi(Jets->at(j).Phi()-METPhi)>(TMath::Pi()-0.4))){
      //std::cout<<"found bad muon jet"<<std::endl;
      result=false;
    }
  }


  //reject events with any jet pt>20, |eta|<2.5 NOT matched to a GenJet (w/in DeltaR<0.3) and chfrac < 0.1
  if(result && runOnSignalMC)
  for(unsigned j = 0; j < Jets->size(); ++j){
    if(Jets->at(j).Pt() <= 20 || fabs(Jets->at(j).Eta())>=2.5) continue;
    bool genMatched = false;
    for(unsigned g = 0; g < GenJets->size(); ++g){
      if(GenJets->at(g).DeltaR(Jets->at(j)) < 0.3) {
         genMatched = true;
         break;
      }
    }
    if(!genMatched && Jets_chargedHadronEnergyFraction->at(j) < 0.1){
      result = false;
      break;
    }
  }

  return result;
}