#define TFMaker_cxx

#include "TFMaker.h"

void TFMaker::Begin(TTree * /*tree*/)
{
    // The Begin() function is called at the start of the query.
    // When running with PROOF Begin() is only called on the client.
    // The tree argument is deprecated (on PROOF 0 is passed).

    TH1::SetDefaultSumw2();
}

void TFMaker::SlaveBegin(TTree * /*tree*/)
{
    // The SlaveBegin() function is called after the Begin() function.
    // When running with PROOF SlaveBegin() is called on each slave server.
    // The tree argument is deprecated (on PROOF 0 is passed).

    if(nicePublication){
        std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
        std::cout<<"!Option for nice plotting enabled! You cannot use Efficiencies.root, as bins might be cut off etc!"<<std::endl;
        std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
    }

    SearchBins_ = new SearchBins(!nicePublication);
    SearchBins_BTags_ = new SearchBins(!nicePublication);

    bTagBins = {0, 0, 0, 0};

    // Initialize Histograms
    TH1::SetDefaultSumw2();
    unsigned nSB = SearchBins_->GetNbins();
    h_CR_SB = new TH1D("h_CR_SB", "h_CR_SB", nSB, 0.5, nSB+0.5);
    h_SR_SB = new TH1D("h_SR_SB", "h_SR_SB", nSB, 0.5, nSB+0.5);
    h_0L1L_SB = new TH1D("h_0L1L_SB", "h_0L1L_SB", nSB, 0.5, nSB+0.5);

    h_CR_SF_SB = new TH1D("h_CR_SF_SB", "h_CR_SF_SB", nSB, 0.5, nSB+0.5);
    h_SR_SF_SB = new TH1D("h_SR_SF_SB", "h_SR_SF_SB", nSB, 0.5, nSB+0.5);
    h_0L1L_SF_SB = new TH1D("h_0L1L_SF_SB", "h_0L1L_SF_SB", nSB, 0.5, nSB+0.5);

    // Use those histograms per sample. You don't want net negative weights
    h_CR_SB_copy = new TH1D("h_CR_SB_copy", "h_CR_SB_copy", nSB, 0.5, nSB+0.5);
    h_SR_SB_copy = new TH1D("h_SR_SB_copy", "h_SR_SB_copy", nSB, 0.5, nSB+0.5);

    h_CR_SF_SB_copy = new TH1D("h_CR_SF_SB_copy", "h_CR_SF_SB_copy", nSB, 0.5, nSB+0.5);
    h_SR_SF_SB_copy = new TH1D("h_SR_SF_SB_copy", "h_SR_SF_SB_copy", nSB, 0.5, nSB+0.5);

    GetOutputList()->Add(h_CR_SB);
    GetOutputList()->Add(h_SR_SB);
    GetOutputList()->Add(h_0L1L_SB);
    GetOutputList()->Add(h_CR_SF_SB);
    GetOutputList()->Add(h_SR_SF_SB);
    GetOutputList()->Add(h_0L1L_SF_SB);

    GetOutputList()->Add(h_CR_SB_copy);
    GetOutputList()->Add(h_SR_SB_copy);
    GetOutputList()->Add(h_CR_SF_SB_copy);
    GetOutputList()->Add(h_SR_SF_SB_copy);

    std::cout<<"----------------"<<std::endl;
    std::cout<<"DeltaPhi Cut: "<<useDeltaPhiCut<<std::endl;
    std::cout<<"----------------"<<std::endl;
}

Bool_t TFMaker::Process(Long64_t entry)
{
    resetValues();

    fChain->GetTree()->GetEntry(entry);

    //  if(entry % 3 != 0) return kTRUE;

    //if(HTgen_cut > 0.01) if(madHT > HTgen_cut) return kTRUE;

    if(HT<minHT_ || MHT< minMHT_ || NJets < minNJets_  ) return kTRUE;
    if(useDeltaPhiCut == 1) if(DeltaPhi1 < deltaPhi1_ || DeltaPhi2 < deltaPhi2_ || DeltaPhi3 < deltaPhi3_ || DeltaPhi4 < deltaPhi4_) return kTRUE;
    if(useDeltaPhiCut == -1) if(!(DeltaPhi1 < deltaPhi1_ || DeltaPhi2 < deltaPhi2_ || DeltaPhi3 < deltaPhi3_ || DeltaPhi4 < deltaPhi4_)) return kTRUE;
    if(applyFilters &&  !FiltersPass() ) return kTRUE;

    GenMuonsNum_ = GenMuons->size();
    GenElectronsNum_ = GenElectrons->size();
    MuonsNum_ = Muons->size();
    ElectronsNum_ = Electrons->size();

    if(GenMuonsNum_ + GenElectronsNum_ == 0) return kTRUE;

    Bin_ = SearchBins_->GetBinNumber(HT,MHT,NJets,BTags);
    if(Bin_ > 900) return kTRUE;

    if(doTopPtReweighting){
        if(GenParticles->size() != GenParticles_PdgId->size()){
            std::cout << "Cannot do top-pT reweighting!"<< std::endl; 
        }else{
            for(unsigned iGen = 0; iGen < GenParticles->size(); iGen++){
                if(std::abs(GenParticles_PdgId->at(iGen)) == 6){
                  topPt.push_back(GenParticles->at(iGen).Pt());
                }
            }

            // https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting#Example
            // Numbers outdated! Use latest numbers from twiki
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
        Weight *= topPtSF;
    }

    TString currentTree = TString(fChain->GetCurrentFile()->GetName());

    if(currentTree != treeName){
        treeName = currentTree;

        TObjArray *optionArray = currentTree.Tokenize("/");
        currFileName = ((TObjString *)(optionArray->At(optionArray->GetEntries()-1)))->String();

        // Make sure you don't have negative number of events per sample
        PushHist(h_CR_SB_copy, h_CR_SB);
        PushHist(h_SR_SB_copy, h_SR_SB);
        PushHist(h_CR_SF_SB_copy, h_CR_SF_SB);
        PushHist(h_SR_SF_SB_copy, h_SR_SF_SB);

        // Open histograms for SFs
        if(SFCR_histFile!=0 || SFSR_histFile!=0){
            h_el_SFCR_etaPt = 0;
            h_el_SFCR_SB = 0;
            h_el_SFSR_etaPt = 0;
            h_el_SFSR_SB = 0;

            h_mu_SFCR_etaPt = 0;
            h_mu_SFCR_SB = 0;
            h_mu_SFSR_etaPt = 0;
            h_mu_SFSR_SB = 0;

            h_di_SFCR_SB = 0;
            h_di_SFSR_SB = 0;

            SFCR_histFile->Close();
            SFCR_histFile = 0;
            SFSR_histFile->Close();
            SFSR_histFile = 0;
        }

        TString SFCR_histFile_path = "";
        TString SFSR_histFile_path = "";
        if((std::string(currFileName.Data()).find(std::string("TTJets"))) != std::string::npos){
            SFCR_histFile_path = "SFCR_0.root";
            SFSR_histFile_path = "SFSR_0.root";
        }else if((std::string(currFileName.Data()).find(std::string("WJetsToLNu"))) != std::string::npos){
            SFCR_histFile_path = "SFCR_1.root";
            SFSR_histFile_path = "SFSR_1.root";
        }else if((std::string(currFileName.Data()).find(std::string("ST_"))) != std::string::npos){
            SFCR_histFile_path = "SFCR_2.root";
            SFSR_histFile_path = "SFSR_2.root";
        }else{
            SFCR_histFile_path = "SFCR_3.root";
            SFSR_histFile_path = "SFSR_3.root";
        }

        SFCR_histFile = TFile::Open(SFCR_histFile_path, "READ");
        SFSR_histFile = TFile::Open(SFSR_histFile_path, "READ");
        h_el_SFCR_etaPt = (TH2D*) SFCR_histFile->Get("h_el_SFCR_etaPt")->Clone();
        h_el_SFCR_SB = (TH1D*) SFCR_histFile->Get("h_el_SFCR_SB")->Clone();
        h_mu_SFCR_etaPt = (TH2D*) SFCR_histFile->Get("h_mu_SFCR_etaPt")->Clone();
        h_mu_SFCR_SB = (TH1D*) SFCR_histFile->Get("h_mu_SFCR_SB")->Clone();
        h_di_SFCR_SB = (TH1D*) SFCR_histFile->Get("h_di_SFCR_SB")->Clone();

        h_el_SFSR_etaPt = (TH2D*) SFSR_histFile->Get("h_el_SFSR_etaPt")->Clone();
        h_el_SFSR_SB = (TH1D*) SFSR_histFile->Get("h_el_SFSR_SB")->Clone();        
        h_mu_SFSR_etaPt = (TH2D*) SFSR_histFile->Get("h_mu_SFSR_etaPt")->Clone();
        h_mu_SFSR_SB = (TH1D*) SFSR_histFile->Get("h_mu_SFSR_SB")->Clone();
        h_di_SFSR_SB = (TH1D*) SFSR_histFile->Get("h_di_SFSR_SB")->Clone();

        if(h_di_SFCR_SB->GetNbinsX() < 100){
            useCombinedBinsCR = true;
            std::cout<<"Using combined bins in CR"<<std::endl;
        }
        if(h_di_SFSR_SB->GetNbinsX() < 100){
            useCombinedBinsSR = true;
            std::cout<<"Using combined bins in SR"<<std::endl;
        }

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

            TFile *skimFile = TFile::Open(path_toSkims+currFileName, "READ");
            btagcorr->SetEffs(skimFile);

            btagcorr->SetCalib(path_bTagCalib);        
            //if(runOnSignalMC){
            //  btagcorr->SetCalibFastSim(path_bTagCalibFastSim);
            //  btagcorr->SetFastSim(true);
            //}
            //else
            btagcorr->SetFastSim(false);
        }

        /*if(runOnSignalMC){
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
        }*/
    }

    /*if(runOnSignalMC){
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
    */
    if(doISRcorr){
        w_isr = isrcorr->GetCorrection(NJetsISR);
        Weight *= w_isr;
    }

    if(doBTagCorr){
    	bTagProb = btagcorr->GetCorrections(Jets,Jets_hadronFlavor,Jets_HTMask);
        bTagBins = {SearchBins_BTags_->GetBinNumber(HT,MHT,NJets,0), SearchBins_BTags_->GetBinNumber(HT,MHT,NJets,1), SearchBins_BTags_->GetBinNumber(HT,MHT,NJets,2), NJets < 3 ? 999 : SearchBins_BTags_->GetBinNumber(HT,MHT,NJets,3)};  
    }
    else{
    	bTagProb = {1, 0, 0, 0};
    	bTagBins = {Bin_, 0, 0, 0};
    }

    //if(runOnSignalMC){
        //  //Account for efficiency of JetID since we cannot apply it on fastSim
        // Weight *= 0.99;
    //}

    //if(useTriggerEffWeight){
        //  if(runOnSignalMC){
        //    Weight *= GetSignalTriggerEffWeight(MHT);
        //  }else{
        //    Weight *= GetTriggerEffWeight(MHT);
        //  }
    //}

    if(doPUreweighting){
        w_pu = puhist->GetBinContent(puhist->GetXaxis()->FindBin(min(TrueNumInteractions,puhist->GetBinLowEdge(puhist->GetNbinsX()+1))));
        Weight *= w_pu;
    }

// We are only interested in the GenLeptons that pass the acceptance, including isotrk veto! (i.e. abs(eta)<2.5, pT>5)
    for(unsigned i=0; i< GenElectrons->size(); i++){
        if(GenElectrons->at(i).Pt() > 5. && std::abs(GenElectrons->at(i).Eta()) < 2.5){
            GenElectronsAcc.push_back(GenElectrons->at(i));
        }
    }
    for(unsigned i=0; i< GenMuons->size(); i++){
        if(GenMuons->at(i).Pt() > 5. && std::abs(GenMuons->at(i).Eta()) < 2.5){
            GenMuonsAcc.push_back(GenMuons->at(i));
        }
    }

    GenMuonsAccNum_ = GenMuonsAcc.size();
    GenElectronsAccNum_ = GenElectronsAcc.size();

    //Define some helpful variables
    if(GenMuonsAccNum_ > 0){
        GenMuonsAccPt_ = GenMuonsAcc.at(0).Pt();
        GenMuonsAccEta_ = GenMuonsAcc.at(0).Eta();
        if(GenMuonsAccNum_ > 1){
            GenMuonsAccPt2_ = GenMuonsAcc.at(1).Pt();
            GenMuonsAccEta2_ = GenMuonsAcc.at(1).Eta();
        }
    }
    if(GenElectronsAccNum_ > 0){
        GenElectronsAccPt_ = GenElectronsAcc.at(0).Pt();
        GenElectronsAccEta_ = GenElectronsAcc.at(0).Eta();
        if(GenElectronsAccNum_ > 1){
            GenElectronsAccPt2_ = GenElectronsAcc.at(1).Pt();
            GenElectronsAccEta2_ = GenElectronsAcc.at(1).Eta();
        }
    }

    // get isoTrack collection from full TAP collection.
    isoTracksNum = isoMuonTracksNum + isoPionTracksNum + isoElectronTracksNum;
    for(unsigned i=0; i< TAPElectronTracks->size(); i++){
        if(TAPElectronTracks_trkiso->at(i) < 0.2 && TAPElectronTracks_mT->at(i) < 100){
            isoElectronTracks.push_back(TAPElectronTracks->at(i));
        }
    }
    if(isoElectronTracks.size() != (unsigned) isoElectronTracksNum){
        std::cout << "WARNING! Number of isoElectronTracks is not correct! Skipping event." << std::endl;
        return kTRUE;
    }

    for(unsigned i=0; i< TAPMuonTracks->size(); i++){
        if(TAPMuonTracks_trkiso->at(i) < 0.2 && TAPMuonTracks_mT->at(i) < 100){
            isoMuonTracks.push_back(TAPMuonTracks->at(i));
        }      
    }
    if(isoMuonTracks.size() != (unsigned) isoMuonTracksNum){
        std::cout << "WARNING! Number of isoMuonTracks is not correct! Skipping event." << std::endl;
        return kTRUE;
    }

    for(unsigned i=0; i< TAPPionTracks->size(); i++){
        if(TAPPionTracks_trkiso->at(i) < 0.1 && TAPPionTracks_mT->at(i) < 100){
            isoPionTracks.push_back(TAPPionTracks->at(i));
        }      
    }
    if(isoPionTracks.size() != (unsigned) isoPionTracksNum){
        std::cout << "WARNING! Number of isoPionTracks is not correct! Skipping event." << std::endl;
        return kTRUE;
    }

    // Match iso leptons/tracks to gen leptons
    // Apply SFs only to prompts
    for(unsigned i=0; i< GenElectronsAccNum_; i++){
        bool matched = false;
        for(unsigned j=0; j< ElectronsNum_; j++){
            if(std::abs(GenElectronsAcc.at(i).Pt() - Electrons->at(j).Pt()) / GenElectronsAcc.at(i).Pt() < 0.1
                && deltaR(GenElectronsAcc.at(i).Eta(), GenElectronsAcc.at(i).Phi(), Electrons->at(j).Eta(), Electrons->at(j).Phi()) < 0.03){
                if(ElectronsPromptNum_==0){
                    ElectronsPromptPt_ = GenElectronsAcc.at(i).Pt();
                    ElectronsPromptEta_ = GenElectronsAcc.at(i).Eta();
                    ElectronsPromptMatch_ = i;
                }else if(ElectronsPromptNum_==1){
                    ElectronsPromptPt2_ = GenElectronsAcc.at(i).Pt();
                    ElectronsPromptEta2_ = GenElectronsAcc.at(i).Eta();
                    ElectronsPromptMatch2_ = i;
                }
                matched = true;
                ElectronsPromptNum_++;
                break;
            }
        }
        if(matched) continue;

        for(int j=0; j< isoElectronTracksNum; j++){
            if(std::abs(GenElectronsAcc.at(i).Pt() - isoElectronTracks.at(j).Pt()) / GenElectronsAcc.at(i).Pt() < 0.1
                && deltaR(GenElectronsAcc.at(i).Eta(), GenElectronsAcc.at(i).Phi(), isoElectronTracks.at(j).Eta(), isoElectronTracks.at(j).Phi()) < 0.03){
                if(ElectronTracksPromptNum_==0){
                    ElectronTracksPromptPt_ = GenElectronsAcc.at(i).Pt();
                    ElectronTracksPromptEta_ = GenElectronsAcc.at(i).Eta();
                    ElectronTracksPromptMatch_ = i;
                }else if(ElectronTracksPromptNum_==1){
                    ElectronTracksPromptPt2_ = GenElectronsAcc.at(i).Pt();
                    ElectronTracksPromptEta2_ = GenElectronsAcc.at(i).Eta();
                    ElectronTracksPromptMatch2_ = i;
                }
                ElectronTracksPromptNum_++;
                break;
            }
        }
    }

    for(unsigned i=0; i< GenMuonsAccNum_; i++){
        bool matched = false;
        for(unsigned j=0; j< MuonsNum_; j++){
            if(std::abs(GenMuonsAcc.at(i).Pt() - Muons->at(j).Pt()) / GenMuonsAcc.at(i).Pt() < 0.1
                && deltaR(GenMuonsAcc.at(i).Eta(), GenMuonsAcc.at(i).Phi(), Muons->at(j).Eta(), Muons->at(j).Phi()) < 0.03){
                if(MuonsPromptNum_==0){
                    MuonsPromptPt_ = GenMuonsAcc.at(i).Pt();
                    MuonsPromptEta_ = GenMuonsAcc.at(i).Eta();
                    MuonsPromptMatch_ = i;
                }else if(MuonsPromptNum_==1){
                    MuonsPromptPt2_ = GenMuonsAcc.at(i).Pt();
                    MuonsPromptEta2_ = GenMuonsAcc.at(i).Eta();
                    MuonsPromptMatch2_ = i;
                }
                matched = true;
                MuonsPromptNum_++;
                break;
            }
        }
        if(matched) continue;

        for(int j=0; j< isoMuonTracksNum; j++){
            if(std::abs(GenMuonsAcc.at(i).Pt() - isoMuonTracks.at(j).Pt()) / GenMuonsAcc.at(i).Pt() < 0.1
                && deltaR(GenMuonsAcc.at(i).Eta(), GenMuonsAcc.at(i).Phi(), isoMuonTracks.at(j).Eta(), isoMuonTracks.at(j).Phi()) < 0.03){
                if(MuonTracksPromptNum_==0){
                    MuonTracksPromptPt_ = GenMuonsAcc.at(i).Pt();
                    MuonTracksPromptEta_ = GenMuonsAcc.at(i).Eta();
                    MuonTracksPromptMatch_ = i;
                }else if(MuonTracksPromptNum_==1){
                    MuonTracksPromptPt2_ = GenMuonsAcc.at(i).Pt();
                    MuonTracksPromptEta2_ = GenMuonsAcc.at(i).Eta();
                    MuonTracksPromptMatch2_ = i;
                }
                MuonTracksPromptNum_++;
                break;
            }
        }
    }

    if(GenMuonsAccNum_ < MuonsPromptNum_ + MuonTracksPromptNum_ || GenElectronsAccNum_ < ElectronsPromptNum_ + ElectronTracksPromptNum_){
        std::cout<<"Mu:"<<GenMuonsAccNum_<<"->"<<MuonsPromptNum_<<"+"<<MuonTracksPromptNum_<<std::endl;
        std::cout<<"El:"<<GenElectronsAccNum_<<"->"<<ElectronsPromptNum_<<"+"<<ElectronTracksPromptNum_<<std::endl;
        std::cout<<"Matching not successful. Skipping event."<<std::endl;
        return kTRUE;
    }



    Weight *= scaleFactorWeight;

    int nLoops = 1;
    if(doBTagCorr) nLoops = (NJets == 2 ? 3 : 4);
    for(int i = 0; i < nLoops; i++){
    	double WeightBtagProb = Weight*bTagProb.at(i);
    	unsigned bTagBin = bTagBins.at(i);

    	// CONTROL REGION
        if(ElectronsNum_ + MuonsNum_ == 1){
            double SF = 1;
            // Bug found bz Aditee
            // double binSF = Bin_;
            // should be (?)
            double binSF = bTagBin;
            if(useCombinedBinsCR) binSF = SearchBins_->GetCombinedBinNumber(HT,MHT,NJets);

            if(ElectronsNum_ == 1){
                mtw = Electrons_MTW->at(0);
                if(GenElectronsAccNum_ == 1 && GenMuonsAccNum_ == 0){
                    SF = GetSF(h_el_SFCR_SB, binSF); 
                }else if(applySFforDilep &&((GenElectronsAccNum_ == 2 && GenMuonsAccNum_ == 0) || (GenElectronsAccNum_ == 1 && GenMuonsAccNum_ == 1))){
                    SF = GetSF(h_di_SFCR_SB, binSF)*GetSF(h_di_SFSR_SB, binSF);
                }else{
                    SF = 1;
                }

                // Don't correct for non-prompts
                if(ElectronsPromptNum_==0) SF = 1;
            }
            if(MuonsNum_ == 1){
                mtw = Muons_MTW->at(0);
                if(GenElectronsAccNum_ == 0 && GenMuonsAccNum_ == 1){
                    SF = GetSF(h_mu_SFCR_SB, binSF); 
                }else if(applySFforDilep && ((GenElectronsAccNum_ == 0 && GenMuonsAccNum_ == 2) || (GenElectronsAccNum_ == 1 && GenMuonsAccNum_ == 1))){
                    SF = GetSF(h_di_SFCR_SB, binSF)*GetSF(h_di_SFSR_SB, binSF);
                }else{
                    SF = 1;
                }

                // Don't correct for non-prompts
                if(MuonsPromptNum_==0) SF = 1;
            }

            if(mtw > 100) return kTRUE;
        
            h_CR_SB_copy->Fill(bTagBin, WeightBtagProb);
            h_CR_SF_SB_copy->Fill(bTagBin, WeightBtagProb*SF);
        }

        // SIGNAL REGION
        if(ElectronsNum_ + MuonsNum_ == 0 && isoTracksNum == 0){
            if(GenElectronsNum_ + GenMuonsNum_ == 0) return kTRUE;

            double SF = 1;
            // Bug found bz Aditee
            // double binSF = Bin_;
            // should be (?)
            double binSF = bTagBin;
            if(useCombinedBinsSR) binSF = SearchBins_->GetCombinedBinNumber(HT,MHT,NJets);

            if(GenElectronsAccNum_ == 1 && GenMuonsAccNum_ == 0){
                SF = GetSF(h_el_SFSR_SB, binSF); 
            }else if(GenElectronsAccNum_ == 0 && GenMuonsAccNum_ == 1){
                SF = GetSF(h_mu_SFSR_SB, binSF); 
            }else if(applySFforDilep && (GenElectronsAccNum_ + GenMuonsAccNum_ == 2)){
                SF = GetSF(h_di_SFSR_SB, binSF)*GetSF(h_di_SFSR_SB, binSF); 
            }else{
                SF = 1;
            }

            h_SR_SB_copy->Fill(bTagBin, WeightBtagProb);
            h_SR_SF_SB_copy->Fill(bTagBin, WeightBtagProb*SF);
        }
        
    }

    return kTRUE;
}

void TFMaker::SlaveTerminate()
{
    // The SlaveTerminate() function is called after all entries or objects
    // have been processed. When running with PROOF SlaveTerminate() is called
    // on each slave server.

    //std::cout<<"--- QCD binning ---"<<std::endl;
    //SearchBinsQCD_->PrintUsed();

    std::cout<<"--- Search bins ---"<<std::endl;
    //SearchBins_->PrintUsedCombined();  
    SearchBins_->PrintUsed();  
}

void TFMaker::Terminate()
{
    // The Terminate() function is the last function to be called during
    // a query. It always runs on the client, it can be used to present
    // the results graphically or save the results to file.

    // Draw Options
    gStyle->SetPaintTextFormat("5.2f");
    gStyle->SetStatW(0.1);
    gStyle->SetStatH(0.1);
    gStyle->SetStatY(202);

    gStyle->SetPalette(87);
    gStyle->SetMarkerSize(1.0);
    gStyle->SetNumberContours(255);

    gStyle->SetTextFont(42);
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetStatFont(42);

    TFile *outPutFile = new TFile(fileName,"RECREATE");

    h_CR_SB = dynamic_cast<TH1D*>(GetOutputList()->FindObject("h_CR_SB"));
    h_CR_SF_SB = dynamic_cast<TH1D*>(GetOutputList()->FindObject("h_CR_SF_SB"));
    h_SR_SB = dynamic_cast<TH1D*>(GetOutputList()->FindObject("h_SR_SB"));
    h_SR_SF_SB = dynamic_cast<TH1D*>(GetOutputList()->FindObject("h_SR_SF_SB"));
    h_0L1L_SB = dynamic_cast<TH1D*>(GetOutputList()->FindObject("h_0L1L_SB"));
    h_0L1L_SF_SB = dynamic_cast<TH1D*>(GetOutputList()->FindObject("h_0L1L_SF_SB"));

    h_CR_SB_copy = dynamic_cast<TH1D*>(GetOutputList()->FindObject("h_CR_SB_copy"));
    h_CR_SF_SB_copy = dynamic_cast<TH1D*>(GetOutputList()->FindObject("h_CR_SF_SB_copy"));
    h_SR_SB_copy = dynamic_cast<TH1D*>(GetOutputList()->FindObject("h_SR_SB_copy"));
    h_SR_SF_SB_copy = dynamic_cast<TH1D*>(GetOutputList()->FindObject("h_SR_SF_SB_copy"));

    PushHist(h_CR_SB_copy, h_CR_SB);
    PushHist(h_SR_SB_copy, h_SR_SB);
    PushHist(h_CR_SF_SB_copy, h_CR_SF_SB);
    PushHist(h_SR_SF_SB_copy, h_SR_SF_SB);


    h_0L1L_SB->Reset();
    h_0L1L_SB->Divide(h_SR_SB, h_CR_SB);

    h_0L1L_SF_SB->Reset();
    h_0L1L_SF_SB->Divide(h_SR_SF_SB, h_CR_SF_SB);

    for(int nX = 1; nX <= h_0L1L_SB->GetXaxis()->GetNbins(); ++nX){
        if(h_0L1L_SB->GetBinContent(nX) < 0) std::cout<<"h_0L1L_SB (Bin "<<nX<<") negative value"<<std::endl;
        if(h_0L1L_SF_SB->GetBinContent(nX) < 0) std::cout<<"h_0L1L_SF_SB (Bin "<<nX<<") negative value"<<std::endl;
    }

    SaveEff(h_CR_SB, outPutFile,"CR_SB;Search region bin number;N_{CR}",false,true);
    SaveEff(h_CR_SF_SB, outPutFile,"CR_SF_SB;Search region bin number;N_{CR}^{SF}",false,true);
    SaveEff(h_SR_SB, outPutFile,"SR_SB;Search region bin number;N_{SR}",false,true);
    SaveEff(h_SR_SF_SB, outPutFile,"SR_SF_SB;Search region bin number;N_{SR}^{SF}",false,true);
    SaveEff(h_0L1L_SB, outPutFile,"0L1L_SB;Search region bin number;TF");
    SaveEff(h_0L1L_SF_SB, outPutFile,"0L1L_SB;Search region bin number;TF^{SF}");

    TH1D* h_TFSB_TF = dynamic_cast<TH1D*>(h_0L1L_SB->Clone("TFSB_TF"));
    h_TFSB_TF->Reset();
    h_TFSB_TF->Divide(h_0L1L_SF_SB, h_0L1L_SB);
    SaveEff(h_TFSB_TF, outPutFile,"TFSB_TF;Search region bin number;TF^{SF}/TF");

    outPutFile->Close();

    cout << "Saved output to " << fileName << endl;

}

void TFMaker::PushHist(TH1* h, TH1* f){
    for(int nX = 1; nX <= h->GetXaxis()->GetNbins(); ++nX){
        if(h->GetBinContent(nX) < 0){
            h->SetBinContent(nX, 0);
            h->SetBinError(nX, 0);
        }
    }

    f->Add(h);
    h->Reset();
}

void TFMaker::SaveEff(TH1* h, TFile* oFile, const char* title, bool xlog, bool ylog)
{
    oFile->cd();

    std::string name = std::string(h->GetName());

    TString TITLE(title);
    while(!TITLE.BeginsWith(";")){
        TITLE.Remove(0,1);
    }
    h->SetTitle(TITLE);
    h->SetMarkerSize(2.0);
    h->UseCurrentStyle();

    h->SetTitleOffset(0.9, "X");
    h->SetTitleOffset(0.8, "Y");
    h->SetTitleOffset(0.7, "Z");
    h->SetTitleSize(0.06, "X");
    h->SetTitleSize(0.06, "Y");
    h->SetTitleSize(0.06, "Z");

    h->SetMarkerStyle(20);
    h->SetMarkerSize(0.8);

    TPaveText *pt = new TPaveText(.15,.96,.35,.95, "NDC");
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetTextFont(42);
    pt->SetTextAlign(31);
    pt->SetTextSize(0.04);
    pt->SetMargin(0.);
    pt->AddText("Simulation");

    TPaveText *pt2 = new TPaveText(.80,.96,.90,.95, "NDC");
    pt2->SetBorderSize(0);
    pt2->SetFillColor(0);
    pt2->SetTextFont(42);
    pt2->SetTextAlign(31);
    pt2->SetTextSize(0.04);
    pt2->SetMargin(0.);
    pt2->AddText("(13 TeV)");

    gROOT->SetBatch(true);    
    TCanvas *c1 = new TCanvas("c1","c1",1);
    c1->SetTopMargin(0.06);
    c1->SetBottomMargin(0.12);
    c1->SetLeftMargin(0.12);

    c1->cd();

    h->GetYaxis()->SetRangeUser(0., 2.);

    if(xlog){
      h->GetXaxis()->SetRangeUser(0.001, h->GetXaxis()->GetBinLowEdge(h->GetNbinsX()+1));
      c1->SetLogx();
    }
    if(ylog){
        h->GetYaxis()->SetRangeUser(0.001, h->GetYaxis()->GetBinLowEdge(h->GetNbinsY()+1));
        c1->SetLogy();
    }

    //h->GetYaxis()->SetRangeUser(0.79, 1.31);
    
    h->Draw("P,E1");
    pt->Draw();
    pt2->Draw();

    //if(name.find(std::string("SFCR")) != std::string::npos || name.find(std::string("SFSR")) != std::string::npos){
        //TObjArray *optionArray = currFileName.Tokenize("_.");
        //TString currTreeName = ((TObjString *)(optionArray->At(0)))->String();
        TString currTreeName("all_noNeg");
        c1->SaveAs("TFs/"+currTreeName+"_"+TString(name)+".pdf");
    //}

    delete c1;
    gROOT->SetBatch(false);

    h->Write();
}

bool TFMaker::FiltersPass()
{
    bool result=true;
    if(useFilterData){
        if(HBHENoiseFilter!=1) result=false;
        if(HBHEIsoNoiseFilter!=1) result=false;
        if(EcalDeadCellTriggerPrimitiveFilter!=1) result=false;    
        if(eeBadScFilter!=1) result=false;

        //if(runOnData){
        //    if(!BadChargedCandidateFilter) result=false;
        //    if(!BadPFMuonFilter) result=false;
        //    if(globalTightHalo2016Filter!=1) result=false;
        //}    
    }
    if(NVtx<=0) result=false;

    // Preliminary filters
    if(PFCaloMETRatio>5) result=false;

    for(unsigned j = 0; j < Jets->size(); j++){
        if(TMath::IsNaN(Jets->at(j).Phi()-METPhi)) result=false;
        if(Jets->at(j).Pt()>200 && Jets_muonEnergyFraction->at(j)>0.5 && (TVector2::Phi_mpi_pi(Jets->at(j).Phi()-METPhi)>(TMath::Pi()-0.4))){
            //std::cout<<"found bad muon jet"<<std::endl;
            result=false;
            break;
        }
    }

    //reject events with any jet pt>20, |eta|<2.5 NOT matched to a GenJet (w/in DeltaR<0.3) and chfrac < 0.1
    /*if (runOnSignalMC)
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
    */
    // Do not apply on fastSim samples!
    //if(!runOnSignalMC) if(!JetID) result=false;
    return result;
}

void TFMaker::resetValues()
{
    GenMuonsAccPt_=0;
    GenMuonsAccEta_=0;
    GenElectronsAccPt_=0;
    GenElectronsAccEta_=0;
    GenMuonsAccPt2_=0;
    GenMuonsAccEta2_=0;
    GenElectronsAccPt2_=0;
    GenElectronsAccEta2_=0;

    MuonsPromptNum_=0;
    ElectronsPromptNum_=0;
    MuonTracksPromptNum_=0;
    ElectronTracksPromptNum_=0;

    MuonsPromptPt_=0;
    MuonsPromptEta_=0;
    ElectronsPromptPt_=0;
    ElectronsPromptEta_=0;

    MuonsPromptPt2_=0;
    MuonsPromptEta2_=0;
    ElectronsPromptPt2_=0;
    ElectronsPromptEta2_=0;

    MuonTracksPromptPt_=0;
    MuonTracksPromptEta_=0;
    ElectronTracksPromptPt_=0;
    ElectronTracksPromptEta_=0;

    MuonTracksPromptPt2_=0;
    MuonTracksPromptEta2_=0;
    ElectronTracksPromptPt2_=0;
    ElectronTracksPromptEta2_=0;

    MuonsPromptMatch_=-1;
    ElectronsPromptMatch_=-1;
    MuonsPromptMatch2_=-1;
    ElectronsPromptMatch2_=-1;

    MuonTracksPromptMatch_=-1;
    ElectronTracksPromptMatch_=-1;
    MuonTracksPromptMatch2_=-1;
    ElectronTracksPromptMatch2_=-1;

    GenElectronsAcc.clear();
    GenMuonsAcc.clear();

    mtw = 0;

    isoElectronTracks.clear();
    isoMuonTracks.clear();
    isoPionTracks.clear();
}