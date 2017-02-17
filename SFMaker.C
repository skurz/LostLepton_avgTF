#define SFMaker_cxx

#include "SFMaker.h"

void SFMaker::Begin(TTree * /*tree*/)
{
    // The Begin() function is called at the start of the query.
    // When running with PROOF Begin() is only called on the client.
    // The tree argument is deprecated (on PROOF 0 is passed).

    TH1::SetDefaultSumw2();
}

void SFMaker::SlaveBegin(TTree * /*tree*/)
{
    // The SlaveBegin() function is called after the Begin() function.
    // When running with PROOF SlaveBegin() is called on each slave server.
    // The tree argument is deprecated (on PROOF 0 is passed).

    SearchBins_ = new SearchBins(true);
    SearchBins_BTags_ = new SearchBins(true);

    bTagBins = {0, 0, 0, 0};

    // Initialize Histograms
    TH1::SetDefaultSumw2();
    unsigned nSB = SearchBins_->GetNbins();
    if(useCombinedBins) nSB = SearchBins_->GetNbinsCombined();
    h_el_nOnePrompt_etaPt = new TH2D("h_el_nOnePrompt_etaPt", "h_el_nOnePrompt_etaPt", nBins_etaElec-1, bins_etaElec, nBins_pT-1, bins_pT);
    h_el_nOnePrompt_SB = new TH1D("h_el_nOnePrompt_SB", "h_el_nOnePrompt_SB", nSB, 0.5, nSB+0.5);    
    h_el_nFoundOnePrompt_etaPt = new TH2D("h_el_nFoundOnePrompt_etaPt", "h_el_nFoundOnePrompt_etaPt", nBins_etaElec-1, bins_etaElec, nBins_pT-1, bins_pT);
    h_el_nFoundOnePrompt_SB = new TH1D("h_el_nFoundOnePrompt_SB", "h_el_nFoundOnePrompt_SB", nSB, 0.5, nSB+0.5);
    h_el_nFoundOnePrompt_SF_etaPt = new TH2D("h_el_nFoundOnePrompt_SF_etaPt", "h_el_nFoundOnePrompt_SF_etaPt", nBins_etaElec-1, bins_etaElec, nBins_pT-1, bins_pT);
    h_el_nFoundOnePrompt_SF_SB = new TH1D("h_el_nFoundOnePrompt_SF_SB", "h_el_nFoundOnePrompt_SF_SB", nSB, 0.5, nSB+0.5);
    h_el_nLostOnePrompt_etaPt = new TH2D("h_el_nLostOnePrompt_etaPt", "h_el_nLostOnePrompt_etaPt", nBins_etaElec-1, bins_etaElec, nBins_pT-1, bins_pT);
    h_el_nLostOnePrompt_SB = new TH1D("h_el_nLostOnePrompt_SB", "h_el_nLostOnePrompt_SB", nSB, 0.5, nSB+0.5);

    h_el_SFCR_etaPt = new TH2D("h_el_SFCR_etaPt", "h_el_SFCR_etaPt", nBins_etaElec-1, bins_etaElec, nBins_pT-1, bins_pT);
    h_el_SFCR_SB = new TH1D("h_el_SFCR_SB", "h_el_SFCR_SB", nSB, 0.5, nSB+0.5);
    h_el_SFSR_etaPt = new TH2D("h_el_SFSR_etaPt", "h_el_SFSR_etaPt", nBins_etaElec-1, bins_etaElec, nBins_pT-1, bins_pT);
    h_el_SFSR_SB = new TH1D("h_el_SFSR_SB", "h_el_SFSR_SB", nSB, 0.5, nSB+0.5);

    h_mu_nOnePrompt_etaPt = new TH2D("h_mu_nOnePrompt_etaPt", "h_mu_nOnePrompt_etaPt", nBins_etaElec-1, bins_etaElec, nBins_pT-1, bins_pT);
    h_mu_nOnePrompt_SB = new TH1D("h_mu_nOnePrompt_SB", "h_mu_nOnePrompt_SB", nSB, 0.5, nSB+0.5);    
    h_mu_nFoundOnePrompt_etaPt = new TH2D("h_mu_nFoundOnePrompt_etaPt", "h_mu_nFoundOnePrompt_etaPt", nBins_etaElec-1, bins_etaElec, nBins_pT-1, bins_pT);
    h_mu_nFoundOnePrompt_SB = new TH1D("h_mu_nFoundOnePrompt_SB", "h_mu_nFoundOnePrompt_SB", nSB, 0.5, nSB+0.5);
    h_mu_nFoundOnePrompt_SF_etaPt = new TH2D("h_mu_nFoundOnePrompt_SF_etaPt", "h_mu_nFoundOnePrompt_SF_etaPt", nBins_etaElec-1, bins_etaElec, nBins_pT-1, bins_pT);
    h_mu_nFoundOnePrompt_SF_SB = new TH1D("h_mu_nFoundOnePrompt_SF_SB", "h_mu_nFoundOnePrompt_SF_SB", nSB, 0.5, nSB+0.5);
    h_mu_nLostOnePrompt_etaPt = new TH2D("h_mu_nLostOnePrompt_etaPt", "h_mu_nLostOnePrompt_etaPt", nBins_etaElec-1, bins_etaElec, nBins_pT-1, bins_pT);
    h_mu_nLostOnePrompt_SB = new TH1D("h_mu_nLostOnePrompt_SB", "h_mu_nLostOnePrompt_SB", nSB, 0.5, nSB+0.5);

    h_mu_SFCR_etaPt = new TH2D("h_mu_SFCR_etaPt", "h_mu_SFCR_etaPt", nBins_etaElec-1, bins_etaElec, nBins_pT-1, bins_pT);
    h_mu_SFCR_SB = new TH1D("h_mu_SFCR_SB", "h_mu_SFCR_SB", nSB, 0.5, nSB+0.5);
    h_mu_SFSR_etaPt = new TH2D("h_mu_SFSR_etaPt", "h_mu_SFSR_etaPt", nBins_etaElec-1, bins_etaElec, nBins_pT-1, bins_pT);
    h_mu_SFSR_SB = new TH1D("h_mu_SFSR_SB", "h_mu_SFSR_SB", nSB, 0.5, nSB+0.5);

    h_di_nTwoPrompt_SB = new TH1D("h_di_nTwoPrompt_SB", "h_di_nTwoPrompt_SB", nSB, 0.5, nSB+0.5);
    h_di_nOneFoundTwoPrompt_SB = new TH1D("h_di_nOneFoundTwoPrompt_SB", "h_di_nOneFoundTwoPrompt_SB", nSB, 0.5, nSB+0.5);
    h_di_nOneFoundTwoPrompt_SF_SB = new TH1D("h_di_nOneFoundTwoPrompt_SF_SB", "h_di_nOneFoundTwoPrompt_SF_SB", nSB, 0.5, nSB+0.5);
    h_di_nTwoFoundTwoPrompt_SB = new TH1D("h_di_nTwoFoundTwoPrompt_SB", "h_di_nTwoFoundTwoPrompt_SB", nSB, 0.5, nSB+0.5);
    h_di_nTwoFoundTwoPrompt_SF_SB = new TH1D("h_di_nTwoFoundTwoPrompt_SF_SB", "h_di_nTwoFoundTwoPrompt_SF_SB", nSB, 0.5, nSB+0.5);
    h_di_nLostTwoPrompt_SB = new TH1D("h_di_nLostTwoPrompt_SB", "h_di_nLostTwoPrompt_SB", nSB, 0.5, nSB+0.5);

    h_di_SFCR_SB = new TH1D("h_di_SFCR_SB", "h_di_SFCR_SB", nSB, 0.5, nSB+0.5);
    h_di_SFSR_SB = new TH1D("h_di_SFSR_SB", "h_di_SFSR_SB", nSB, 0.5, nSB+0.5);

    GetOutputList()->Add(h_el_nOnePrompt_etaPt);
    GetOutputList()->Add(h_el_nOnePrompt_SB);
    GetOutputList()->Add(h_el_nFoundOnePrompt_etaPt);
    GetOutputList()->Add(h_el_nFoundOnePrompt_SB);
    GetOutputList()->Add(h_el_nFoundOnePrompt_SF_etaPt);
    GetOutputList()->Add(h_el_nFoundOnePrompt_SF_SB);
    GetOutputList()->Add(h_el_nLostOnePrompt_etaPt);
    GetOutputList()->Add(h_el_nLostOnePrompt_SB);
    GetOutputList()->Add(h_el_SFCR_etaPt);
    GetOutputList()->Add(h_el_SFCR_SB);
    GetOutputList()->Add(h_el_SFSR_etaPt);
    GetOutputList()->Add(h_el_SFSR_SB);

    GetOutputList()->Add(h_mu_nOnePrompt_etaPt);
    GetOutputList()->Add(h_mu_nOnePrompt_SB);
    GetOutputList()->Add(h_mu_nFoundOnePrompt_etaPt);
    GetOutputList()->Add(h_mu_nFoundOnePrompt_SB);
    GetOutputList()->Add(h_mu_nFoundOnePrompt_SF_etaPt);
    GetOutputList()->Add(h_mu_nFoundOnePrompt_SF_SB);
    GetOutputList()->Add(h_mu_nLostOnePrompt_etaPt);
    GetOutputList()->Add(h_mu_nLostOnePrompt_SB);
    GetOutputList()->Add(h_mu_SFCR_etaPt);
    GetOutputList()->Add(h_mu_SFCR_SB);
    GetOutputList()->Add(h_mu_SFSR_etaPt);
    GetOutputList()->Add(h_mu_SFSR_SB);

    GetOutputList()->Add(h_di_nTwoPrompt_SB);
    GetOutputList()->Add(h_di_nOneFoundTwoPrompt_SB);
    GetOutputList()->Add(h_di_nOneFoundTwoPrompt_SF_SB);
    GetOutputList()->Add(h_di_nTwoFoundTwoPrompt_SB);
    GetOutputList()->Add(h_di_nTwoFoundTwoPrompt_SF_SB);
    GetOutputList()->Add(h_di_nLostTwoPrompt_SB);
    GetOutputList()->Add(h_di_SFCR_SB);
    GetOutputList()->Add(h_di_SFSR_SB);

    std::cout<<"----------------"<<std::endl;
    std::cout<<"DeltaPhi Cut: "<<useDeltaPhiCut<<std::endl;
    std::cout<<"----------------"<<std::endl;
}

Bool_t SFMaker::Process(Long64_t entry)
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

    if(GenMuonsNum_ + GenElectronsNum_ == 0) return kTRUE;

    if(useCombinedBins){
        Bin_ = SearchBins_->GetCombinedBinNumber(HT,MHT,NJets);
    }else{
        Bin_ = SearchBins_->GetBinNumber(HT,MHT,NJets,BTags);
    }    
    if(Bin_ > 900) return kTRUE;

    // TH1 cannot properly deal with negative bin contents
    // At most 1% difference in SFs expected (rare BGs only)
    //if(Weight < 0) return kTRUE;

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
        if(useCombinedBins){
            bTagBins = {Bin_, Bin_, Bin_, Bin_};
        }else{
            bTagBins = {SearchBins_BTags_->GetBinNumber(HT,MHT,NJets,0), SearchBins_BTags_->GetBinNumber(HT,MHT,NJets,1), SearchBins_BTags_->GetBinNumber(HT,MHT,NJets,2), NJets < 3 ? 999 : SearchBins_BTags_->GetBinNumber(HT,MHT,NJets,3)};
        }  
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
    if(GenMuonsAccNum_ + GenElectronsAccNum_ == 0) return kTRUE;

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

    ElectronsNum_ = Electrons->size();
    MuonsNum_ = Muons->size();

    // get isoTrack collection from full TAP collection
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

    int nLoops = 1;
    if(doBTagCorr) nLoops = (NJets == 2 ? 3 : 4);
    for(int i = 0; i < nLoops; i++){
    	double WeightBtagProb = Weight*bTagProb.at(i);
    	unsigned bTagBin = bTagBins.at(i);

    	if(GenMuonsAccNum_ == 1 && GenElectronsAccNum_ == 0){
    		h_mu_nOnePrompt_etaPt->Fill(GenMuonsAccEta_, GenMuonsAccPt_, WeightBtagProb);
	    	h_mu_nOnePrompt_SB->Fill(bTagBin, WeightBtagProb);

            isoSF = GetSF(h_muIsoSF, GenMuonsAccPt_, std::abs(GenMuonsAccEta_));
            recoSF = GetSF(h_muIDSF, GenMuonsAccPt_, std::abs(GenMuonsAccEta_));
            if(GenMuonsAccPt_ > 10) trackingSF = GetSF(h_muTrkHighPtSF, GenMuonsAccEta_);
            else trackingSF = GetSF(h_muTrkLowPtSF, GenMuonsAccEta_);

    		if(MuonsPromptNum_ == 1){
                double WeightCorr = WeightBtagProb * isoSF * recoSF * trackingSF;

		        h_mu_nFoundOnePrompt_etaPt->Fill(GenMuonsAccEta_, GenMuonsAccPt_, WeightBtagProb);
			    h_mu_nFoundOnePrompt_SB->Fill(bTagBin, WeightBtagProb);
			    h_mu_nFoundOnePrompt_SF_etaPt->Fill(GenMuonsAccEta_, GenMuonsAccPt_, WeightCorr);
			    h_mu_nFoundOnePrompt_SF_SB->Fill(bTagBin, WeightCorr);
    		}else if(includeIsotrkVeto && MuonsPromptNum_ == 0 && MuonTracksPromptNum_ == 1){
                double WeightCorr = WeightBtagProb * trackingSF;

                h_mu_nFoundOnePrompt_etaPt->Fill(GenMuonsAccEta_, GenMuonsAccPt_, WeightBtagProb);
                h_mu_nFoundOnePrompt_SB->Fill(bTagBin, WeightBtagProb);
                h_mu_nFoundOnePrompt_SF_etaPt->Fill(GenMuonsAccEta_, GenMuonsAccPt_, WeightCorr);
                h_mu_nFoundOnePrompt_SF_SB->Fill(bTagBin, WeightCorr);
            }else if(MuonsPromptNum_ == 0 && (!includeIsotrkVeto || MuonTracksPromptNum_ == 0)){
    			h_mu_nLostOnePrompt_etaPt->Fill(GenMuonsAccEta_, GenMuonsAccPt_, WeightBtagProb);
	    		h_mu_nLostOnePrompt_SB->Fill(bTagBin, WeightBtagProb);
    		}else{
                std::cout<<"SingleMu: "<<MuonsPromptNum_<<"+"<<MuonTracksPromptNum_<<std::endl;
            }
    	}

    	if(GenMuonsAccNum_ == 0 && GenElectronsAccNum_ == 1){
    		h_el_nOnePrompt_etaPt->Fill(GenElectronsAccEta_, GenElectronsAccPt_, WeightBtagProb);
	    	h_el_nOnePrompt_SB->Fill(bTagBin, WeightBtagProb);

            isoSF = GetSF(h_elecIsoSF, GenElectronsAccPt_, std::abs(GenElectronsAccEta_));
            recoSF = GetSF(h_elecIDSF, GenElectronsAccPt_, std::abs(GenElectronsAccEta_));
            trackingSF = GetSF(h_elecTrkSF, GenElectronsAccEta_, GenElectronsAccPt_); 

    		if(ElectronsPromptNum_ == 1){
                double WeightCorr = WeightBtagProb * isoSF * recoSF * trackingSF;

                h_el_nFoundOnePrompt_etaPt->Fill(GenElectronsAccEta_, GenElectronsAccPt_, WeightBtagProb);
			    h_el_nFoundOnePrompt_SB->Fill(bTagBin, WeightBtagProb);
			    h_el_nFoundOnePrompt_SF_etaPt->Fill(GenElectronsAccEta_, GenElectronsAccPt_, WeightCorr);
			    h_el_nFoundOnePrompt_SF_SB->Fill(bTagBin, WeightCorr);
    		}else if(includeIsotrkVeto && ElectronsPromptNum_ == 0 && ElectronTracksPromptNum_ == 1){
                double WeightCorr = WeightBtagProb * trackingSF;

                h_el_nFoundOnePrompt_etaPt->Fill(GenElectronsAccEta_, GenElectronsAccPt_, WeightBtagProb);
                h_el_nFoundOnePrompt_SB->Fill(bTagBin, WeightBtagProb);
                h_el_nFoundOnePrompt_SF_etaPt->Fill(GenElectronsAccEta_, GenElectronsAccPt_, WeightCorr);
                h_el_nFoundOnePrompt_SF_SB->Fill(bTagBin, WeightCorr);
            }else if(ElectronsPromptNum_ == 0 && (!includeIsotrkVeto || ElectronTracksPromptNum_ == 0)){
    			h_el_nLostOnePrompt_etaPt->Fill(GenElectronsAccEta_, GenElectronsAccPt_, WeightBtagProb);
	    		h_el_nLostOnePrompt_SB->Fill(bTagBin, WeightBtagProb);
    		}else{
                std::cout<<"SingleElec: "<<ElectronsPromptNum_<<"+"<<ElectronTracksPromptNum_<<std::endl;
            }	
    	}

        if(GenMuonsAccNum_ == 2 && GenElectronsAccNum_ == 0){
            h_di_nTwoPrompt_SB->Fill(bTagBin, WeightBtagProb);

            isoSF = GetSF(h_muIsoSF, GenMuonsAccPt_, std::abs(GenMuonsAccEta_));
            recoSF = GetSF(h_muIDSF, GenMuonsAccPt_, std::abs(GenMuonsAccEta_));
            if(GenMuonsAccPt_ > 10) trackingSF = GetSF(h_muTrkHighPtSF, GenMuonsAccEta_);
            else trackingSF = GetSF(h_muTrkLowPtSF, GenMuonsAccEta_);

            isoSF2 = GetSF(h_muIsoSF, GenMuonsAccPt2_, std::abs(GenMuonsAccEta2_));
            recoSF2 = GetSF(h_muIDSF, GenMuonsAccPt2_, std::abs(GenMuonsAccEta2_));
            if(GenMuonsAccPt2_ > 10) trackingSF = GetSF(h_muTrkHighPtSF, GenMuonsAccEta2_);
            else trackingSF2 = GetSF(h_muTrkLowPtSF, GenMuonsAccEta2_);

            if(MuonsPromptNum_ == 2){
                double WeightCorr = WeightBtagProb * isoSF * recoSF * trackingSF * isoSF2 * recoSF2 * trackingSF2;

                h_di_nTwoFoundTwoPrompt_SB->Fill(bTagBin, WeightBtagProb);
                h_di_nTwoFoundTwoPrompt_SF_SB->Fill(bTagBin, WeightCorr);
            }else if(includeIsotrkVeto && MuonsPromptNum_ == 1 && MuonTracksPromptNum_ == 1){
                double WeightCorr = 1;
                if(MuonsPromptMatch_ == 0)  WeightCorr = WeightBtagProb * isoSF * recoSF * trackingSF * trackingSF2;
                else WeightCorr = WeightBtagProb * trackingSF * isoSF2 * recoSF2 * trackingSF2;

                h_di_nTwoFoundTwoPrompt_SB->Fill(bTagBin, WeightBtagProb);
                h_di_nTwoFoundTwoPrompt_SF_SB->Fill(bTagBin, WeightCorr);
            }else if(includeIsotrkVeto && MuonsPromptNum_ == 0 && MuonTracksPromptNum_ == 2){
                double WeightCorr = WeightBtagProb * trackingSF * trackingSF2;

                h_di_nTwoFoundTwoPrompt_SB->Fill(bTagBin, WeightBtagProb);
                h_di_nTwoFoundTwoPrompt_SF_SB->Fill(bTagBin, WeightCorr);
            }else if(MuonsPromptNum_ == 1 && (!includeIsotrkVeto || MuonTracksPromptNum_ == 0)){
                double WeightCorr = 1;
                if(MuonsPromptMatch_ == 0)  WeightCorr = WeightBtagProb * isoSF * recoSF * trackingSF;
                else WeightCorr = WeightBtagProb * trackingSF * isoSF2 * recoSF2 * trackingSF2;

                h_di_nOneFoundTwoPrompt_SB->Fill(bTagBin, WeightBtagProb);
                h_di_nOneFoundTwoPrompt_SF_SB->Fill(bTagBin, WeightCorr);
            }else if(includeIsotrkVeto && MuonsPromptNum_ == 0 && MuonTracksPromptNum_ == 1){
                double WeightCorr = 1;
                if(MuonsPromptMatch_ == 0)  WeightCorr = WeightBtagProb * trackingSF;
                else WeightCorr = WeightBtagProb * trackingSF2;

                h_di_nOneFoundTwoPrompt_SB->Fill(bTagBin, WeightBtagProb);
                h_di_nOneFoundTwoPrompt_SF_SB->Fill(bTagBin, WeightCorr);
            }else if(MuonsPromptNum_ == 0 && (!includeIsotrkVeto || MuonTracksPromptNum_ == 0)){
                h_di_nLostTwoPrompt_SB->Fill(bTagBin, WeightBtagProb);
            }else{
                std::cout<<"DiMu: "<<GenMuonsAccNum_<<"->"<<MuonsPromptNum_<<"+"<<MuonTracksPromptNum_<<std::endl;
            }       
        }

        if(GenMuonsAccNum_ == 0 && GenElectronsAccNum_ == 2){
            h_di_nTwoPrompt_SB->Fill(bTagBin, WeightBtagProb);

            isoSF = GetSF(h_elecIsoSF, GenElectronsAccPt_, std::abs(GenElectronsAccEta_));
            recoSF = GetSF(h_elecIDSF, GenElectronsAccPt_, std::abs(GenElectronsAccEta_));
            trackingSF = GetSF(h_elecTrkSF, GenElectronsAccEta_, GenElectronsAccPt_); 

            isoSF2 = GetSF(h_elecIsoSF, GenElectronsAccPt2_, std::abs(GenElectronsAccEta2_));
            recoSF2 = GetSF(h_elecIDSF, GenElectronsAccPt2_, std::abs(GenElectronsAccEta2_));
            trackingSF2 = GetSF(h_elecTrkSF, GenElectronsAccEta2_, GenElectronsAccPt2_); 

            if(ElectronsPromptNum_ == 2){
                double WeightCorr = WeightBtagProb * isoSF * recoSF * trackingSF * isoSF2 * recoSF2 * trackingSF2;

                h_di_nTwoFoundTwoPrompt_SB->Fill(bTagBin, WeightBtagProb);
                h_di_nTwoFoundTwoPrompt_SF_SB->Fill(bTagBin, WeightCorr);
            }else if(includeIsotrkVeto && ElectronsPromptNum_ == 1 && ElectronTracksPromptNum_ == 1){
                double WeightCorr = 1;
                if(MuonsPromptMatch_ == 0)  WeightCorr = WeightBtagProb * isoSF * recoSF * trackingSF * trackingSF2;
                else WeightCorr = WeightBtagProb * trackingSF * isoSF2 * recoSF2 * trackingSF2;

                h_di_nTwoFoundTwoPrompt_SB->Fill(bTagBin, WeightBtagProb);
                h_di_nTwoFoundTwoPrompt_SF_SB->Fill(bTagBin, WeightCorr);
            }else if(includeIsotrkVeto && ElectronsPromptNum_ == 0 && ElectronTracksPromptNum_ == 2){
                double WeightCorr = WeightBtagProb * trackingSF * trackingSF2;

                h_di_nTwoFoundTwoPrompt_SB->Fill(bTagBin, WeightBtagProb);
                h_di_nTwoFoundTwoPrompt_SF_SB->Fill(bTagBin, WeightCorr);
            }else if(ElectronsPromptNum_ == 1 && (!includeIsotrkVeto || ElectronTracksPromptNum_ == 0)){
                double WeightCorr = 1;
                if(MuonsPromptMatch_ == 0)  WeightCorr = WeightBtagProb * isoSF * recoSF * trackingSF;
                else WeightCorr = WeightBtagProb * trackingSF * isoSF2 * recoSF2 * trackingSF2;

                h_di_nOneFoundTwoPrompt_SB->Fill(bTagBin, WeightBtagProb);
                h_di_nOneFoundTwoPrompt_SF_SB->Fill(bTagBin, WeightCorr);
            }else if(includeIsotrkVeto && ElectronsPromptNum_ == 0 && ElectronTracksPromptNum_ == 1){
                double WeightCorr = 1;
                if(MuonsPromptMatch_ == 0)  WeightCorr = WeightBtagProb * trackingSF;
                else WeightCorr = WeightBtagProb * trackingSF2;

                h_di_nOneFoundTwoPrompt_SB->Fill(bTagBin, WeightBtagProb);
                h_di_nOneFoundTwoPrompt_SF_SB->Fill(bTagBin, WeightCorr);
            }else if(ElectronsPromptNum_ == 0 && (!includeIsotrkVeto || ElectronTracksPromptNum_ == 0)){
                h_di_nLostTwoPrompt_SB->Fill(bTagBin, WeightBtagProb);
            }else{
                std::cout<<"DiElec: "<<GenElectronsAccNum_<<"->"<<ElectronsPromptNum_<<"+"<<ElectronTracksPromptNum_<<std::endl;
            }   
        }

        if(GenMuonsAccNum_ == 1 && GenElectronsAccNum_ == 1){
            h_di_nTwoPrompt_SB->Fill(bTagBin, WeightBtagProb);

            isoSF = GetSF(h_muIsoSF, GenMuonsAccPt_, std::abs(GenMuonsAccEta_));
            recoSF = GetSF(h_muIDSF, GenMuonsAccPt_, std::abs(GenMuonsAccEta_));
            if(GenMuonsAccPt_ > 10) trackingSF = GetSF(h_muTrkHighPtSF, GenMuonsAccEta_);
            else trackingSF = GetSF(h_muTrkLowPtSF, GenMuonsAccEta_);

            isoSF2 = GetSF(h_elecIsoSF, GenElectronsAccPt_, std::abs(GenElectronsAccEta_));
            recoSF2 = GetSF(h_elecIDSF, GenElectronsAccPt_, std::abs(GenElectronsAccEta_));
            trackingSF2 = GetSF(h_elecTrkSF, GenElectronsAccEta_, GenElectronsAccPt_); 

            if(ElectronsPromptNum_ == 1 && MuonsPromptNum_ == 1){
                double WeightCorr = WeightBtagProb * isoSF * recoSF * trackingSF * isoSF2 * recoSF2 * trackingSF2;

                h_di_nTwoFoundTwoPrompt_SB->Fill(bTagBin, WeightBtagProb);
                h_di_nTwoFoundTwoPrompt_SF_SB->Fill(bTagBin, WeightCorr);
            }else if(includeIsotrkVeto && ((MuonsPromptNum_ == 1 && ElectronTracksPromptNum_ == 1) || (ElectronsPromptNum_ == 1 && MuonTracksPromptNum_ == 1))){
                double WeightCorr = WeightBtagProb * isoSF * recoSF * trackingSF * trackingSF2;
                if(ElectronsPromptNum_ == 1) WeightCorr = WeightBtagProb * trackingSF * isoSF2 * recoSF2 * trackingSF2;

                h_di_nTwoFoundTwoPrompt_SB->Fill(bTagBin, WeightBtagProb);
                h_di_nTwoFoundTwoPrompt_SF_SB->Fill(bTagBin, WeightCorr);
            }else if(includeIsotrkVeto && MuonsPromptNum_ == 0 && ElectronsPromptNum_ == 0 && MuonTracksPromptNum_ == 1 && ElectronTracksPromptNum_ == 1){
                double WeightCorr = WeightBtagProb * trackingSF * trackingSF2;

                h_di_nTwoFoundTwoPrompt_SB->Fill(bTagBin, WeightBtagProb);
                h_di_nTwoFoundTwoPrompt_SF_SB->Fill(bTagBin, WeightCorr);
            }else if((MuonsPromptNum_ == 1 && ElectronsPromptNum_ == 0 && (!includeIsotrkVeto || ElectronTracksPromptNum_ == 0)) || (ElectronsPromptNum_ == 1 && MuonsPromptNum_ == 0 && (!includeIsotrkVeto || MuonTracksPromptNum_ == 0))){
                double WeightCorr = WeightBtagProb * isoSF * recoSF * trackingSF;
                if(ElectronsPromptNum_ == 1) WeightCorr = WeightBtagProb * isoSF2 * recoSF2 * trackingSF2;

                h_di_nOneFoundTwoPrompt_SB->Fill(bTagBin, WeightBtagProb);
                h_di_nOneFoundTwoPrompt_SF_SB->Fill(bTagBin, WeightCorr);
            }else if(includeIsotrkVeto && ((MuonsPromptNum_ == 0 && ElectronsPromptNum_ == 0 && ElectronTracksPromptNum_ == 0 && MuonTracksPromptNum_ == 1) || (MuonsPromptNum_ == 0 && ElectronsPromptNum_ == 0 && MuonTracksPromptNum_ == 0 && ElectronTracksPromptNum_ == 1))){
                double WeightCorr = WeightBtagProb * trackingSF;
                if(ElectronTracksPromptNum_ == 1) WeightCorr = WeightBtagProb * trackingSF2;

                h_di_nOneFoundTwoPrompt_SB->Fill(bTagBin, WeightBtagProb);
                h_di_nOneFoundTwoPrompt_SF_SB->Fill(bTagBin, WeightCorr);
            }else if(MuonsPromptNum_ == 0 && ElectronsPromptNum_ == 0 && (!includeIsotrkVeto || (MuonTracksPromptNum_ == 0 && ElectronTracksPromptNum_ == 0))){
                h_di_nLostTwoPrompt_SB->Fill(bTagBin, WeightBtagProb);
            }else{
                std::cout<<"DiMuEl: "<<GenMuonsAccNum_<<"/"<<GenElectronsAccNum_<<"->"<<MuonsPromptNum_<<"+"<<MuonTracksPromptNum_<<"/"<<ElectronsPromptNum_<<"+"<<ElectronTracksPromptNum_<<std::endl;
            }
        }
    }

    return kTRUE;
}

void SFMaker::SlaveTerminate()
{
    // The SlaveTerminate() function is called after all entries or objects
    // have been processed. When running with PROOF SlaveTerminate() is called
    // on each slave server.

    //std::cout<<"--- QCD binning ---"<<std::endl;
    //SearchBinsQCD_->PrintUsed();

    std::cout<<"--- Search bins ---"<<std::endl;
    if(useCombinedBins){
        SearchBins_->PrintUsedCombined();  
    }else{
        SearchBins_->PrintUsed();  
    }
}

void SFMaker::Terminate()
{
    // The Terminate() function is the last function to be called during
    // a query. It always runs on the client, it can be used to present
    // the results graphically or save the results to file.

    // Draw Options
    gStyle->SetPaintTextFormat("5.2f");
    gStyle->SetStatW(0.1);
    gStyle->SetStatH(0.1);
    gStyle->SetStatY(202);
    gStyle->SetTitleYOffset(1.3);

    gStyle->SetPalette(87);
    gStyle->SetMarkerSize(1.3);

    TFile *outPutFile = new TFile(fileName,"RECREATE");

    h_el_nOnePrompt_etaPt = dynamic_cast<TH2D*>(GetOutputList()->FindObject("h_el_nOnePrompt_etaPt"));
    h_el_nOnePrompt_SB = dynamic_cast<TH1D*>(GetOutputList()->FindObject("h_el_nOnePrompt_SB"));
    h_el_nFoundOnePrompt_etaPt = dynamic_cast<TH2D*>(GetOutputList()->FindObject("h_el_nFoundOnePrompt_etaPt"));
    h_el_nFoundOnePrompt_SB = dynamic_cast<TH1D*>(GetOutputList()->FindObject("h_el_nFoundOnePrompt_SB"));
    h_el_nFoundOnePrompt_SF_etaPt = dynamic_cast<TH2D*>(GetOutputList()->FindObject("h_el_nFoundOnePrompt_SF_etaPt"));
    h_el_nFoundOnePrompt_SF_SB = dynamic_cast<TH1D*>(GetOutputList()->FindObject("h_el_nFoundOnePrompt_SF_SB"));
    h_el_nLostOnePrompt_etaPt = dynamic_cast<TH2D*>(GetOutputList()->FindObject("h_el_nLostOnePrompt_etaPt"));
    h_el_nLostOnePrompt_SB = dynamic_cast<TH1D*>(GetOutputList()->FindObject("h_el_nLostOnePrompt_SB"));

    h_el_SFCR_etaPt = dynamic_cast<TH2D*>(GetOutputList()->FindObject("h_el_SFCR_etaPt"));
    h_el_SFCR_SB = dynamic_cast<TH1D*>(GetOutputList()->FindObject("h_el_SFCR_SB"));

    h_mu_nOnePrompt_etaPt = dynamic_cast<TH2D*>(GetOutputList()->FindObject("h_mu_nOnePrompt_etaPt"));
    h_mu_nOnePrompt_SB = dynamic_cast<TH1D*>(GetOutputList()->FindObject("h_mu_nOnePrompt_SB"));
    h_mu_nFoundOnePrompt_etaPt = dynamic_cast<TH2D*>(GetOutputList()->FindObject("h_mu_nFoundOnePrompt_etaPt"));
    h_mu_nFoundOnePrompt_SB = dynamic_cast<TH1D*>(GetOutputList()->FindObject("h_mu_nFoundOnePrompt_SB"));
    h_mu_nFoundOnePrompt_SF_etaPt = dynamic_cast<TH2D*>(GetOutputList()->FindObject("h_mu_nFoundOnePrompt_SF_etaPt"));
    h_mu_nFoundOnePrompt_SF_SB = dynamic_cast<TH1D*>(GetOutputList()->FindObject("h_mu_nFoundOnePrompt_SF_SB"));
    h_mu_nLostOnePrompt_etaPt = dynamic_cast<TH2D*>(GetOutputList()->FindObject("h_mu_nLostOnePrompt_etaPt"));
    h_mu_nLostOnePrompt_SB = dynamic_cast<TH1D*>(GetOutputList()->FindObject("h_mu_nLostOnePrompt_SB"));

    h_mu_SFCR_etaPt = dynamic_cast<TH2D*>(GetOutputList()->FindObject("h_mu_SFCR_etaPt"));
    h_mu_SFCR_SB = dynamic_cast<TH1D*>(GetOutputList()->FindObject("h_mu_SFCR_SB"));

    h_el_SFSR_etaPt = dynamic_cast<TH2D*>(GetOutputList()->FindObject("h_el_SFSR_etaPt"));
    h_el_SFSR_SB = dynamic_cast<TH1D*>(GetOutputList()->FindObject("h_el_SFSR_SB"));
    h_mu_SFSR_etaPt = dynamic_cast<TH2D*>(GetOutputList()->FindObject("h_mu_SFSR_etaPt"));
    h_mu_SFSR_SB = dynamic_cast<TH1D*>(GetOutputList()->FindObject("h_mu_SFSR_SB"));

    h_di_nTwoPrompt_SB = dynamic_cast<TH1D*>(GetOutputList()->FindObject("h_di_nTwoPrompt_SB"));
	h_di_nOneFoundTwoPrompt_SB = dynamic_cast<TH1D*>(GetOutputList()->FindObject("h_di_nOneFoundTwoPrompt_SB"));
	h_di_nOneFoundTwoPrompt_SF_SB = dynamic_cast<TH1D*>(GetOutputList()->FindObject("h_di_nOneFoundTwoPrompt_SF_SB"));
	h_di_nTwoFoundTwoPrompt_SB = dynamic_cast<TH1D*>(GetOutputList()->FindObject("h_di_nTwoFoundTwoPrompt_SB"));
	h_di_nTwoFoundTwoPrompt_SF_SB = dynamic_cast<TH1D*>(GetOutputList()->FindObject("h_di_nTwoFoundTwoPrompt_SF_SB"));
	h_di_nLostTwoPrompt_SB = dynamic_cast<TH1D*>(GetOutputList()->FindObject("h_di_nLostTwoPrompt_SB"));

	h_di_SFCR_SB = dynamic_cast<TH1D*>(GetOutputList()->FindObject("h_di_SFCR_SB"));
	h_di_SFSR_SB = dynamic_cast<TH1D*>(GetOutputList()->FindObject("h_di_SFSR_SB"));

    for(int nX = 1; nX <= h_el_nOnePrompt_SB->GetXaxis()->GetNbins(); ++nX){
        if(h_el_nOnePrompt_SB->GetBinContent(nX) < 0){
            h_el_nOnePrompt_SB->SetBinContent(nX, 0);
            std::cout<<"h_el_nOnePrompt_SB (Bin "<<nX<<") negative value"<<std::endl;
        }
        if(h_el_nFoundOnePrompt_SB->GetBinContent(nX) < 0){
            h_el_nFoundOnePrompt_SB->SetBinContent(nX, 0);
            std::cout<<"h_el_nFoundOnePrompt_SB (Bin "<<nX<<") negative value"<<std::endl;
        }
        if(h_el_nFoundOnePrompt_SF_SB->GetBinContent(nX) < 0){
            h_el_nFoundOnePrompt_SF_SB->SetBinContent(nX, 0);
            std::cout<<"h_el_nFoundOnePrompt_SF_SB (Bin "<<nX<<") negative value"<<std::endl;
        }
        if(h_el_nLostOnePrompt_SB->GetBinContent(nX) < 0){
            h_el_nLostOnePrompt_SB->SetBinContent(nX, 0);
            std::cout<<"h_el_nLostOnePrompt_SB (Bin "<<nX<<") negative value"<<std::endl;
        }
        if(h_mu_nOnePrompt_SB->GetBinContent(nX) < 0){
            h_mu_nOnePrompt_SB->SetBinContent(nX, 0);
            std::cout<<"h_mu_nOnePrompt_SB (Bin "<<nX<<") negative value"<<std::endl;
        }
        if(h_mu_nFoundOnePrompt_SB->GetBinContent(nX) < 0){
            h_mu_nFoundOnePrompt_SB->SetBinContent(nX, 0);
            std::cout<<"h_mu_nFoundOnePrompt_SB (Bin "<<nX<<") negative value"<<std::endl;
        }
        if(h_mu_nFoundOnePrompt_SF_SB->GetBinContent(nX) < 0){
            h_mu_nFoundOnePrompt_SF_SB->SetBinContent(nX, 0);
            std::cout<<"h_mu_nFoundOnePrompt_SF_SB (Bin "<<nX<<") negative value"<<std::endl;
        }
        if(h_mu_nLostOnePrompt_SB->GetBinContent(nX) < 0){
            h_mu_nLostOnePrompt_SB->SetBinContent(nX, 0);
            std::cout<<"h_mu_nLostOnePrompt_SB (Bin "<<nX<<") negative value"<<std::endl;
        }        
        if(h_di_nTwoPrompt_SB->GetBinContent(nX) < 0){
            h_di_nTwoPrompt_SB->SetBinContent(nX, 0);
            std::cout<<"h_di_nTwoPrompt_SB (Bin "<<nX<<") negative value"<<std::endl;
        }
        if(h_di_nOneFoundTwoPrompt_SB->GetBinContent(nX) < 0){
            h_di_nOneFoundTwoPrompt_SB->SetBinContent(nX, 0);
            std::cout<<"h_di_nOneFoundTwoPrompt_SB (Bin "<<nX<<") negative value"<<std::endl;
        }
        if(h_di_nOneFoundTwoPrompt_SF_SB->GetBinContent(nX) < 0){
            h_di_nOneFoundTwoPrompt_SF_SB->SetBinContent(nX, 0);
            std::cout<<"h_di_nOneFoundTwoPrompt_SF_SB (Bin "<<nX<<") negative value"<<std::endl;
        }
        if(h_di_nTwoFoundTwoPrompt_SB->GetBinContent(nX) < 0){
            h_di_nTwoFoundTwoPrompt_SB->SetBinContent(nX, 0);
            std::cout<<"h_di_nTwoFoundTwoPrompt_SB (Bin "<<nX<<") negative value"<<std::endl;
        }
        if(h_di_nTwoFoundTwoPrompt_SF_SB->GetBinContent(nX) < 0){
            h_di_nTwoFoundTwoPrompt_SF_SB->SetBinContent(nX, 0);
            std::cout<<"h_di_nTwoFoundTwoPrompt_SF_SB (Bin "<<nX<<") negative value"<<std::endl;
        }
        if(h_di_nLostTwoPrompt_SB->GetBinContent(nX) < 0){
            h_di_nLostTwoPrompt_SB->SetBinContent(nX, 0);
            std::cout<<"h_di_nLostTwoPrompt_SB (Bin "<<nX<<") negative value"<<std::endl;
        }
    }



    ///////////////
    /// Singleleptonic
    ////////////

    // Without SFs:
    // nOnePrompt = nFoundOnePrompt + nLostOnePrompt
    // nOnePrompt is constant (only depends on xsec)

    // With SFs:
    // nOnePrompt = SFCR * nFoundOnePrompt + SFSR * nLostOnePrompt

    // Solve for SFSR:
    // SFSR = (nOnePrompt - SFCR * nFoundOnePrompt) / nLostOnePrompt

    h_el_SFSR_etaPt->Reset();
    h_el_SFSR_etaPt->Add(h_el_nOnePrompt_etaPt);
    h_el_SFSR_etaPt->Add(h_el_nFoundOnePrompt_SF_etaPt, -1);
    h_el_SFSR_etaPt->Divide(h_el_nLostOnePrompt_etaPt);

    h_el_SFSR_SB->Reset();
    h_el_SFSR_SB->Add(h_el_nOnePrompt_SB);
    h_el_SFSR_SB->Add(h_el_nFoundOnePrompt_SF_SB, -1);
    h_el_SFSR_SB->Divide(h_el_nLostOnePrompt_SB);

    h_mu_SFSR_etaPt->Reset();
    h_mu_SFSR_etaPt->Add(h_mu_nOnePrompt_etaPt);
    h_mu_SFSR_etaPt->Add(h_mu_nFoundOnePrompt_SF_etaPt, -1);
    h_mu_SFSR_etaPt->Divide(h_mu_nLostOnePrompt_etaPt);

    h_mu_SFSR_SB->Reset();
    h_mu_SFSR_SB->Add(h_mu_nOnePrompt_SB);
    h_mu_SFSR_SB->Add(h_mu_nFoundOnePrompt_SF_SB, -1);
    h_mu_SFSR_SB->Divide(h_mu_nLostOnePrompt_SB);


    // SF for CRs
    h_mu_SFCR_etaPt->Divide(h_mu_nFoundOnePrompt_SF_etaPt, h_mu_nFoundOnePrompt_etaPt);
	h_mu_SFCR_SB->Divide(h_mu_nFoundOnePrompt_SF_SB, h_mu_nFoundOnePrompt_SB);

	h_el_SFCR_etaPt->Divide(h_el_nFoundOnePrompt_SF_etaPt, h_el_nFoundOnePrompt_etaPt);
	h_el_SFCR_SB->Divide(h_el_nFoundOnePrompt_SF_SB, h_el_nFoundOnePrompt_SB);


    // Save histograms
    //SaveEff(h_el_nOnePrompt_etaPt, outPutFile, false, true);
    //SaveEff(h_el_nOnePrompt_SB, outPutFile);
    h_el_nFoundOnePrompt_etaPt->Divide(h_el_nOnePrompt_etaPt);
    SaveEff(h_el_nFoundOnePrompt_etaPt, outPutFile, false, true);
    h_el_nFoundOnePrompt_SB->Divide(h_el_nOnePrompt_SB);
    SaveEff(h_el_nFoundOnePrompt_SB, outPutFile);
    h_el_nFoundOnePrompt_SF_etaPt->Divide(h_el_nOnePrompt_etaPt);
    SaveEff(h_el_nFoundOnePrompt_SF_etaPt, outPutFile, false, true);
    h_el_nFoundOnePrompt_SF_SB->Divide(h_el_nOnePrompt_SB);
    SaveEff(h_el_nFoundOnePrompt_SF_SB, outPutFile);
    h_el_nLostOnePrompt_etaPt->Divide(h_el_nOnePrompt_etaPt);
    SaveEff(h_el_nLostOnePrompt_etaPt, outPutFile, false, true);
    h_el_nLostOnePrompt_SB->Divide(h_el_nOnePrompt_SB);
    SaveEff(h_el_nLostOnePrompt_SB, outPutFile);

    for(int nX = 1; nX <= h_el_SFCR_SB->GetXaxis()->GetNbins(); ++nX){
        h_el_SFCR_SB->SetBinError(nX, 0);
        h_el_SFSR_SB->SetBinError(nX, 0);

        if(h_el_SFCR_SB->GetBinContent(nX) < 1e-8) h_el_SFCR_SB->SetBinContent(nX, 1);
        if(h_el_SFSR_SB->GetBinContent(nX) < 1e-8) h_el_SFSR_SB->SetBinContent(nX, 1);

        // Fix for sample with negative weights
        if(h_el_SFCR_SB->GetBinContent(nX) > 1) h_el_SFCR_SB->SetBinContent(nX, 1);
        if(h_el_SFSR_SB->GetBinContent(nX) < 1) h_el_SFSR_SB->SetBinContent(nX, 1);
    }

    SaveEff(h_el_SFCR_etaPt, outPutFile, false, true);
    SaveEff(h_el_SFCR_SB, outPutFile);

    SaveEff(h_el_SFSR_etaPt, outPutFile, false, true);
    SaveEff(h_el_SFSR_SB, outPutFile);

    //SaveEff(h_mu_nOnePrompt_etaPt, outPutFile, false, true);
    //SaveEff(h_mu_nOnePrompt_SB, outPutFile);
    h_mu_nFoundOnePrompt_etaPt->Divide(h_mu_nOnePrompt_etaPt);
    SaveEff(h_mu_nFoundOnePrompt_etaPt, outPutFile, false, true);
    h_mu_nFoundOnePrompt_SB->Divide(h_mu_nOnePrompt_SB);
    SaveEff(h_mu_nFoundOnePrompt_SB, outPutFile);
    h_mu_nFoundOnePrompt_SF_etaPt->Divide(h_mu_nOnePrompt_etaPt);
    SaveEff(h_mu_nFoundOnePrompt_SF_etaPt, outPutFile, false, true);
    h_mu_nFoundOnePrompt_SF_SB->Divide(h_mu_nOnePrompt_SB);
    SaveEff(h_mu_nFoundOnePrompt_SF_SB, outPutFile);
    h_mu_nLostOnePrompt_etaPt->Divide(h_mu_nOnePrompt_etaPt);
    SaveEff(h_mu_nLostOnePrompt_etaPt, outPutFile, false, true);
    h_mu_nLostOnePrompt_SB->Divide(h_mu_nOnePrompt_SB);
    SaveEff(h_mu_nLostOnePrompt_SB, outPutFile);

    for(int nX = 1; nX <= h_mu_SFCR_SB->GetXaxis()->GetNbins(); ++nX){
        h_mu_SFCR_SB->SetBinError(nX, 0);
        h_mu_SFSR_SB->SetBinError(nX, 0);
        
        if(h_mu_SFCR_SB->GetBinContent(nX) < 1e-8) h_mu_SFCR_SB->SetBinContent(nX, 1);
        if(h_mu_SFSR_SB->GetBinContent(nX) < 1e-8) h_mu_SFSR_SB->SetBinContent(nX, 1);

        // Fix for sample with negative weights
        if(h_mu_SFCR_SB->GetBinContent(nX) > 1) h_mu_SFCR_SB->SetBinContent(nX, 1);
        if(h_mu_SFSR_SB->GetBinContent(nX) < 1) h_mu_SFSR_SB->SetBinContent(nX, 1);
    }

    SaveEff(h_mu_SFCR_etaPt, outPutFile, false, true);
    SaveEff(h_mu_SFCR_SB, outPutFile);

    SaveEff(h_mu_SFSR_etaPt, outPutFile, false, true);
    SaveEff(h_mu_SFSR_SB, outPutFile);


    ///////////////
    /// Dileptonic
    ////////////

    // Similar as singleleptonic case
    // With SFs:
    // nTwoPrompt = SFCR^2 * nTwoFoundTwoPrompt + SFCR * SFSR * nOneFoundTwoPrompt + SFSR^2 * nLostTwoPrompt
    // => 0 = SFSR^2 * nLostTwoPrompt + SFSR * (SFCR * nOneFoundTwoPrompt) + ((SFCR^2 * nTwoFoundTwoPrompt) - nTwoPrompt)

    // Solve for SFSR: (quadratic equation with only 1 positive solution)
    // SFSR = (-nOneFoundTwoPrompt + sqrt((nOneFoundTwoPrompt)^2 - 4 * nLostTwoPrompt * (nTwoFoundTwoPrompt - nTwoPrompt))) /  (2 * nLostTwoPrompt)
    // Better: Use Muller's method (smaller rounding errors)
    // x = 2c / (-b - sqrt(b^2 - 4ac))
    // SFSR = 2 * (nTwoFoundTwoPrompt - nTwoPrompt) / (-nOneFoundTwoPrompt - sqrt(nOneFoundTwoPrompt^2 - 4 * nLostTwoPrompt * (nTwoFoundTwoPrompt - nTwoPrompt)))

	for(int nX = 1; nX <= h_di_SFSR_SB->GetXaxis()->GetNbins(); ++nX){
		double SFSR_dilep = 1;
        if(h_di_nTwoPrompt_SB->GetBinContent(nX) > 0 && h_di_nTwoFoundTwoPrompt_SF_SB->GetBinContent(nX) - h_di_nTwoPrompt_SB->GetBinContent(nX) < 0 && (h_di_nLostTwoPrompt_SB->GetBinContent(nX) > 0 || h_di_nOneFoundTwoPrompt_SF_SB->GetBinContent(nX) > 0))
            SFSR_dilep = (2. * (h_di_nTwoFoundTwoPrompt_SF_SB->GetBinContent(nX) - h_di_nTwoPrompt_SB->GetBinContent(nX))) / (-h_di_nOneFoundTwoPrompt_SF_SB->GetBinContent(nX) - std::sqrt(h_di_nOneFoundTwoPrompt_SF_SB->GetBinContent(nX)*h_di_nOneFoundTwoPrompt_SF_SB->GetBinContent(nX) - 4. * h_di_nLostTwoPrompt_SB->GetBinContent(nX) * (h_di_nTwoFoundTwoPrompt_SF_SB->GetBinContent(nX) - h_di_nTwoPrompt_SB->GetBinContent(nX))));
		//else if(h_di_nLostTwoPrompt_SB->GetBinContent(nX) > 0)
        //SFSR_dilep = (-h_di_nOneFoundTwoPrompt_SF_SB->GetBinContent(nX) + std::sqrt(h_di_nOneFoundTwoPrompt_SF_SB->GetBinContent(nX)*h_di_nOneFoundTwoPrompt_SF_SB->GetBinContent(nX) - 4. * h_di_nLostTwoPrompt_SB->GetBinContent(nX) * (h_di_nTwoFoundTwoPrompt_SF_SB->GetBinContent(nX) - h_di_nTwoPrompt_SB->GetBinContent(nX)))) / (2. * h_di_nLostTwoPrompt_SB->GetBinContent(nX));
		h_di_SFSR_SB->SetBinContent(nX, SFSR_dilep);

        //std::cout<<SFSR_dilep<<"   "<<h_di_nTwoPrompt_SB->GetBinContent(nX)<<"  "<<h_di_nTwoFoundTwoPrompt_SF_SB->GetBinContent(nX)<<"/"<<h_di_nTwoFoundTwoPrompt_SB->GetBinContent(nX)<<"  "<<h_di_nOneFoundTwoPrompt_SF_SB->GetBinContent(nX)<<"/"<<h_di_nOneFoundTwoPrompt_SB->GetBinContent(nX)<<"    "<<h_di_nLostTwoPrompt_SB->GetBinContent(nX)<<std::endl;

		double SFCR_dilep = 1;
		if(h_di_nTwoFoundTwoPrompt_SB->GetBinContent(nX) > 0) SFCR_dilep = std::sqrt(h_di_nTwoFoundTwoPrompt_SF_SB->GetBinContent(nX)/h_di_nTwoFoundTwoPrompt_SB->GetBinContent(nX));
		h_di_SFCR_SB->SetBinContent(nX, SFCR_dilep);

        // Fix for sample with negative weights
        if(h_di_SFCR_SB->GetBinContent(nX) > 1) h_di_SFCR_SB->SetBinContent(nX, 1);
        if(h_di_SFSR_SB->GetBinContent(nX) < 1) h_di_SFSR_SB->SetBinContent(nX, 1);
	}

    h_di_nTwoFoundTwoPrompt_SB->Divide(h_di_nTwoPrompt_SB);
    SaveEff(h_di_nTwoFoundTwoPrompt_SB, outPutFile);
    h_di_nOneFoundTwoPrompt_SB->Divide(h_di_nTwoPrompt_SB);
    SaveEff(h_di_nOneFoundTwoPrompt_SB, outPutFile);
    h_di_nTwoFoundTwoPrompt_SF_SB->Divide(h_di_nTwoPrompt_SB);
    SaveEff(h_di_nTwoFoundTwoPrompt_SF_SB, outPutFile);
    h_di_nOneFoundTwoPrompt_SF_SB->Divide(h_di_nTwoPrompt_SB);
    SaveEff(h_di_nOneFoundTwoPrompt_SF_SB, outPutFile);
    h_di_nLostTwoPrompt_SB->Divide(h_di_nTwoPrompt_SB);
    SaveEff(h_di_nLostTwoPrompt_SB, outPutFile);

    SaveEff(h_di_SFCR_SB, outPutFile);
    SaveEff(h_di_SFSR_SB, outPutFile);

    outPutFile->Close();

    cout << "Saved output to " << fileName << endl;

}

void SFMaker::SaveEff(TH1* h, TFile* oFile, bool xlog, bool ylog)
{
    oFile->cd();

    //h->SetTitle(TString("Simulation, L=3 fb^{-1}, #sqrt{s}=13 TeV ") + TString(title));
    h->SetMarkerSize(2.0);
    h->UseCurrentStyle();

    gROOT->SetBatch(true);    
    TCanvas *c1 = new TCanvas("c1","c1",1);
    c1->cd();
    if(xlog){
      c1->SetLogx();
      h->GetXaxis()->SetRangeUser(0.001, h->GetXaxis()->GetBinLowEdge(h->GetNbinsX()+1));
    }
    if(ylog){
    	c1->SetLogy();
    	h->GetYaxis()->SetRangeUser(0.001, h->GetYaxis()->GetBinLowEdge(h->GetNbinsY()+1));
   	}

    std::string name = std::string(h->GetName());
    if(name.find(std::string("SFCR")) != std::string::npos || name.find(std::string("SFSR")) != std::string::npos){
    	if(name.find(std::string("SB")) != std::string::npos) h->GetYaxis()->SetRangeUser(0.79, 1.31);
    	if(name.find(std::string("etaPt")) != std::string::npos){
            h->GetYaxis()->SetRangeUser(5., 500.);
            h->GetZaxis()->SetRangeUser(0.79, 1.31);
        }
    }else{
    	if(name.find(std::string("SB")) != std::string::npos) h->GetYaxis()->SetRangeUser(0.01, 1.01);
    	if(name.find(std::string("etaPt")) != std::string::npos){
            h->GetYaxis()->SetRangeUser(5., 500.);
            h->GetZaxis()->SetRangeUser(0.01, 1.01);
        }
    }
    
    h->Draw("ColZ,Text");

    //if(name.find(std::string("SFCR")) != std::string::npos || name.find(std::string("SFSR")) != std::string::npos){
    	TObjArray *optionArray = currFileName.Tokenize("_.");
    	TString currTreeName = ((TObjString *)(optionArray->At(0)))->String();
    	c1->SaveAs("SFs/"+currTreeName+"_"+TString(name)+".pdf");
	//}

    delete c1;
    gROOT->SetBatch(false);

    h->Write();
}

bool SFMaker::FiltersPass()
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

void SFMaker::resetValues()
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