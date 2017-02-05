#ifndef SEARCHBINS_H
#define SEARCHBINS_H

#include "TVector2.h"
#include <iostream>
#include <vector>
#include <utility>
#include <TROOT.h>


class Bin
{
public:
	Bin(){}
	Bin(int HTIdx, double HTmin, double HTmax, int MHTIdx, double MHTmin, double MHTmax, int NJetsIdx, int NJetsmin, int NJetsmax, int BTagsIdx, int BTagsmin, int BTagsmax)
	{
		HTIdx_=HTIdx;
		HTmin_=HTmin;
		HTmax_=HTmax;
		MHTIdx_=MHTIdx;
		MHTmin_=MHTmin;
		MHTmax_=MHTmax;
		NJetsIdx_=NJetsIdx;
		NJetsmin_=NJetsmin;
		NJetsmax_=NJetsmax;
		BTagsIdx_=BTagsIdx;
		BTagsmin_=BTagsmin;
		BTagsmax_=BTagsmax;
	}
	~Bin(){}

	double getHTmin(){return HTmin_;}
	double getHTmax(){return HTmax_;}
	double getMHTmin(){return MHTmin_;}
	double getMHTmax(){return MHTmax_;}
	int getNJetsmin(){return NJetsmin_;}
	int getNJetsmax(){return NJetsmax_;}
	int getBTagsmin(){return BTagsmin_;}
	int getBTagsmax(){return BTagsmax_;}

	// return string according to assigned bin numbering scheme (e.g. NJets0_BTags1_HT2_MHT3)
	TString getString();
	std::vector<int> getIdx();
private:
	double HTmin_, HTmax_, MHTmin_, MHTmax_;
	int NJetsmin_, NJetsmax_, BTagsmin_, BTagsmax_;
	// Indixes in the bin numbering scheme. QCD control bins have negative values!
	int NJetsIdx_, BTagsIdx_, HTIdx_, MHTIdx_;
};

TString Bin::getString(){
	if(MHTIdx_ >= 0) return TString::Format("NJets%d_BTags%d_MHT%d_HT%d", NJetsIdx_, BTagsIdx_, MHTIdx_, HTIdx_);
	else return TString::Format("NJets%d_BTags%d_MHTC_HTC%d", NJetsIdx_, BTagsIdx_, (-HTIdx_));
}

std::vector<int> Bin::getIdx(){
	std::vector<int> v;
	v.push_back(NJetsIdx_);
	v.push_back(BTagsIdx_);
	v.push_back(MHTIdx_);
	v.push_back(HTIdx_);
	return v;
}

class SearchBins
{
public:
  SearchBins() : SearchBins(false) {}
  SearchBins(bool useQCDbinning);
  unsigned int GetBinNumber(double HT, double MHT, int NJets, int BTags);
  unsigned int GetCombinedBinNumber(double HT, double MHT, int NJets);
  unsigned int GetNbins();
  unsigned int GetNbinsCombined();
  void PrintUsed();
  void PrintUsedCombined();
  Bin* GetSearchBin(int i);
  
  ~SearchBins(){}
private:
  std::vector<Bin> bins_;
  std::vector<Bin> binsCombined_;
  std::vector<int> usedBin_;
  std::vector<int> usedCombinedBin_;

  std::vector<std::pair<std::pair<double, double>, int>> HT_bins_NJets26;
  std::vector<std::pair<std::pair<double, double>, int>> MHT_bins_NJets26;
  std::vector<std::pair<std::pair<double, double>, int>> HT_bins_NJets7Inf;
  std::vector<std::pair<std::pair<double, double>, int>> MHT_bins_NJets7Inf;
  std::vector<std::pair<std::pair<int, int>, int>> NJets_bins;
  std::vector<std::pair<std::pair<int, int>, int>> BTags_bins_NJets3Inf;
  std::vector<std::pair<std::pair<int, int>, int>> BTags_bins_NJets2;

  Long64_t nSkipped = 0;
};


SearchBins::SearchBins(bool useQCDbinning)
{
	// Add QCD control bins
	// Assign negatice indixes for QCD HTMHT bins
	if(useQCDbinning){
		// [HT range]
		HT_bins_NJets26.push_back(std::make_pair(std::make_pair(300, 500), -1));
		HT_bins_NJets26.push_back(std::make_pair(std::make_pair(500, 1000), -2));
		HT_bins_NJets26.push_back(std::make_pair(std::make_pair(1000, 99999), -3));

		// [MHT range]
		MHT_bins_NJets26.push_back(std::make_pair(std::make_pair(250, 300), -1));
		MHT_bins_NJets26.push_back(std::make_pair(std::make_pair(250, 300), -1));
		MHT_bins_NJets26.push_back(std::make_pair(std::make_pair(250, 300), -1));

		// [HT range]
		HT_bins_NJets7Inf.push_back(std::make_pair(std::make_pair(500, 1000), -2));
		HT_bins_NJets7Inf.push_back(std::make_pair(std::make_pair(1000, 99999), -3));

		// [MHT range]
		MHT_bins_NJets7Inf.push_back(std::make_pair(std::make_pair(250, 300), -1));
		MHT_bins_NJets7Inf.push_back(std::make_pair(std::make_pair(250, 300), -1));
	}

	// Standard Search Bin Definition
	// 10 HTMHT bins
	// [HT range]
	// Search Bins
	HT_bins_NJets26.push_back(std::make_pair(std::make_pair(300, 500), 0));
	HT_bins_NJets26.push_back(std::make_pair(std::make_pair(500, 1000), 1));
	HT_bins_NJets26.push_back(std::make_pair(std::make_pair(1000, 99999), 2));
	HT_bins_NJets26.push_back(std::make_pair(std::make_pair(350, 500), 3));
	HT_bins_NJets26.push_back(std::make_pair(std::make_pair(500, 1000), 4));
	HT_bins_NJets26.push_back(std::make_pair(std::make_pair(1000, 99999), 5));
	HT_bins_NJets26.push_back(std::make_pair(std::make_pair(500, 1000), 6));
	HT_bins_NJets26.push_back(std::make_pair(std::make_pair(1000, 99999), 7));
	HT_bins_NJets26.push_back(std::make_pair(std::make_pair(750, 1500), 8));
	HT_bins_NJets26.push_back(std::make_pair(std::make_pair(1500, 99999), 9));

	// [MHT range]
	MHT_bins_NJets26.push_back(std::make_pair(std::make_pair(300, 350), 0));
	MHT_bins_NJets26.push_back(std::make_pair(std::make_pair(300, 350), 0));
	MHT_bins_NJets26.push_back(std::make_pair(std::make_pair(300, 350), 0));
	MHT_bins_NJets26.push_back(std::make_pair(std::make_pair(350, 500), 1));
	MHT_bins_NJets26.push_back(std::make_pair(std::make_pair(350, 500), 1));
	MHT_bins_NJets26.push_back(std::make_pair(std::make_pair(350, 500), 1));
	MHT_bins_NJets26.push_back(std::make_pair(std::make_pair(500, 750), 2));
	MHT_bins_NJets26.push_back(std::make_pair(std::make_pair(500, 750), 2));
	MHT_bins_NJets26.push_back(std::make_pair(std::make_pair(750, 99999), 3));
	MHT_bins_NJets26.push_back(std::make_pair(std::make_pair(750, 99999) ,3));

	// [HT range]
	HT_bins_NJets7Inf.push_back(std::make_pair(std::make_pair(500, 1000), 1));
	HT_bins_NJets7Inf.push_back(std::make_pair(std::make_pair(1000, 99999), 2));
	HT_bins_NJets7Inf.push_back(std::make_pair(std::make_pair(500, 1000), 4));
	HT_bins_NJets7Inf.push_back(std::make_pair(std::make_pair(1000, 99999), 5));
	HT_bins_NJets7Inf.push_back(std::make_pair(std::make_pair(500, 1000), 6));
	HT_bins_NJets7Inf.push_back(std::make_pair(std::make_pair(1000, 99999), 7));
	HT_bins_NJets7Inf.push_back(std::make_pair(std::make_pair(750, 1500), 8));
	HT_bins_NJets7Inf.push_back(std::make_pair(std::make_pair(1500, 99999), 9));

	// [MHT range]
	MHT_bins_NJets7Inf.push_back(std::make_pair(std::make_pair(300, 350), 0));
	MHT_bins_NJets7Inf.push_back(std::make_pair(std::make_pair(300, 350), 0));
	MHT_bins_NJets7Inf.push_back(std::make_pair(std::make_pair(350, 500), 1));
	MHT_bins_NJets7Inf.push_back(std::make_pair(std::make_pair(350, 500), 1));
	MHT_bins_NJets7Inf.push_back(std::make_pair(std::make_pair(500, 750), 2));
	MHT_bins_NJets7Inf.push_back(std::make_pair(std::make_pair(500, 750), 2));
	MHT_bins_NJets7Inf.push_back(std::make_pair(std::make_pair(750, 99999), 3));
	MHT_bins_NJets7Inf.push_back(std::make_pair(std::make_pair(750, 99999), 3));

	// NJets
	NJets_bins.push_back(std::make_pair(std::make_pair(2, 2), 0));
	NJets_bins.push_back(std::make_pair(std::make_pair(3, 4), 1));
	NJets_bins.push_back(std::make_pair(std::make_pair(5, 6), 2));
	NJets_bins.push_back(std::make_pair(std::make_pair(7, 8), 3));
	NJets_bins.push_back(std::make_pair(std::make_pair(9, 999), 4));

	// BTags (different for NJets==2)
	BTags_bins_NJets3Inf.push_back(std::make_pair(std::make_pair(0, 0), 0));
	BTags_bins_NJets3Inf.push_back(std::make_pair(std::make_pair(1, 1), 1));
	BTags_bins_NJets3Inf.push_back(std::make_pair(std::make_pair(2, 2), 2));
	BTags_bins_NJets3Inf.push_back(std::make_pair(std::make_pair(3, 99), 3));

	BTags_bins_NJets2.push_back(std::make_pair(std::make_pair(0, 0), 0));
	BTags_bins_NJets2.push_back(std::make_pair(std::make_pair(1, 1), 1));
	BTags_bins_NJets2.push_back(std::make_pair(std::make_pair(2, 99), 2));
	

	if(HT_bins_NJets26.size() != MHT_bins_NJets26.size())	std::cout << "SearchBins.h: No valid definition of HTMHT bins!" << std::endl;
	if(HT_bins_NJets7Inf.size() != MHT_bins_NJets7Inf.size())	std::cout << "SearchBins.h: No valid definition of HTMHT bins!" << std::endl;

	// Create bins and push back
	for(unsigned i_NJets = 0; i_NJets < NJets_bins.size(); ++i_NJets){
		// NJets=2 (skip highest BTags bin)
		if(NJets_bins[i_NJets].first.second == 2){
			for(unsigned i_bTags = 0; i_bTags < BTags_bins_NJets2.size(); ++i_bTags){
				for(unsigned i_HTMHT = 0; i_HTMHT < HT_bins_NJets26.size(); ++i_HTMHT){
					bins_.push_back(Bin(HT_bins_NJets26[i_HTMHT].second,	HT_bins_NJets26[i_HTMHT].first.first,	HT_bins_NJets26[i_HTMHT].first.second,
										MHT_bins_NJets26[i_HTMHT].second,	MHT_bins_NJets26[i_HTMHT].first.first,	MHT_bins_NJets26[i_HTMHT].first.second,
										NJets_bins[i_NJets].second,			NJets_bins[i_NJets].first.first,		NJets_bins[i_NJets].first.second,
										BTags_bins_NJets2[i_bTags].second,	BTags_bins_NJets2[i_bTags].first.first, BTags_bins_NJets2[i_bTags].first.second));
				}
			}
		// NJets=3-6 (standard)
		}else if(NJets_bins[i_NJets].first.second <= 6){
			for(unsigned i_bTags = 0; i_bTags < BTags_bins_NJets3Inf.size(); ++i_bTags){
				for(unsigned i_HTMHT = 0; i_HTMHT < HT_bins_NJets26.size(); ++i_HTMHT){
					bins_.push_back(Bin(HT_bins_NJets26[i_HTMHT].second,			HT_bins_NJets26[i_HTMHT].first.first,		HT_bins_NJets26[i_HTMHT].first.second,
										MHT_bins_NJets26[i_HTMHT].second,			MHT_bins_NJets26[i_HTMHT].first.first,		MHT_bins_NJets26[i_HTMHT].first.second,
										NJets_bins[i_NJets].second,					NJets_bins[i_NJets].first.first,			NJets_bins[i_NJets].first.second,
										BTags_bins_NJets3Inf[i_bTags].second, 		BTags_bins_NJets3Inf[i_bTags].first.first,	BTags_bins_NJets3Inf[i_bTags].first.second));
				}
			}
		// NJets=7-Inf (skip lowest HT bins)
		}else{
			for(unsigned i_bTags = 0; i_bTags < BTags_bins_NJets3Inf.size(); ++i_bTags){
				for(unsigned i_HTMHT = 0; i_HTMHT < HT_bins_NJets7Inf.size(); ++i_HTMHT){
					bins_.push_back(Bin(HT_bins_NJets7Inf[i_HTMHT].second,			HT_bins_NJets7Inf[i_HTMHT].first.first,		HT_bins_NJets7Inf[i_HTMHT].first.second,
										MHT_bins_NJets7Inf[i_HTMHT].second,			MHT_bins_NJets7Inf[i_HTMHT].first.first,	MHT_bins_NJets7Inf[i_HTMHT].first.second,
										NJets_bins[i_NJets].second,					NJets_bins[i_NJets].first.first,			NJets_bins[i_NJets].first.second,
										BTags_bins_NJets3Inf[i_bTags].second,		BTags_bins_NJets3Inf[i_bTags].first.first,	BTags_bins_NJets3Inf[i_bTags].first.second));
				}
			}
		}

		///////////////
		// Combine bins in BTag multiplicity
		///////////////
		// NJets=2 (skip highest BTags bin)
		if(NJets_bins[i_NJets].first.second == 2){
			for(unsigned i_HTMHT = 0; i_HTMHT < HT_bins_NJets26.size(); ++i_HTMHT){
				binsCombined_.push_back(Bin(HT_bins_NJets26[i_HTMHT].second,	HT_bins_NJets26[i_HTMHT].first.first,	HT_bins_NJets26[i_HTMHT].first.second,
									MHT_bins_NJets26[i_HTMHT].second,	MHT_bins_NJets26[i_HTMHT].first.first,	MHT_bins_NJets26[i_HTMHT].first.second,
									NJets_bins[i_NJets].second,			NJets_bins[i_NJets].first.first,		NJets_bins[i_NJets].first.second,
									999,								0, 										99));
			}
		// NJets=3-6 (standard)
		}else if(NJets_bins[i_NJets].first.second <= 6){
			for(unsigned i_HTMHT = 0; i_HTMHT < HT_bins_NJets26.size(); ++i_HTMHT){
				binsCombined_.push_back(Bin(HT_bins_NJets26[i_HTMHT].second,			HT_bins_NJets26[i_HTMHT].first.first,		HT_bins_NJets26[i_HTMHT].first.second,
									MHT_bins_NJets26[i_HTMHT].second,			MHT_bins_NJets26[i_HTMHT].first.first,		MHT_bins_NJets26[i_HTMHT].first.second,
									NJets_bins[i_NJets].second,					NJets_bins[i_NJets].first.first,			NJets_bins[i_NJets].first.second,
									999,								0, 										99));
			}
		// NJets=7-Inf (skip lowest HT bins)
		}else{
			for(unsigned i_HTMHT = 0; i_HTMHT < HT_bins_NJets7Inf.size(); ++i_HTMHT){
				binsCombined_.push_back(Bin(HT_bins_NJets7Inf[i_HTMHT].second,			HT_bins_NJets7Inf[i_HTMHT].first.first,		HT_bins_NJets7Inf[i_HTMHT].first.second,
									MHT_bins_NJets7Inf[i_HTMHT].second,			MHT_bins_NJets7Inf[i_HTMHT].first.first,	MHT_bins_NJets7Inf[i_HTMHT].first.second,
									NJets_bins[i_NJets].second,					NJets_bins[i_NJets].first.first,			NJets_bins[i_NJets].first.second,
									999,								0, 										99));
			}
		}

	}

	std::cout<<"Loaded bins: "<<bins_.size()<<"("<<binsCombined_.size()<<")"<<std::endl;
	for(unsigned int i=0; i<bins_.size();i++)
	{
		//std::cout<<"Bin["<<i<<"]: HT["<<bins_[i].HTmin_<<","<<bins_[i].HTmax_<<"] MHT["<<bins_[i].MHTmin_<<","<<bins_[i].MHTmax_<<"] NJets["<<bins_[i].NJetsmin_<<","<<bins_[i].NJetsmax_<<"] BTags["<<bins_[i].BTagsmin_<<","<<bins_[i].BTagsmax_<<"]\n";
	    usedBin_.push_back(0); 
	}
	for(unsigned int i=0; i<binsCombined_.size();i++)
	{
		usedCombinedBin_.push_back(0); 
	}
}

Bin* SearchBins::GetSearchBin(int i){
	if(i>(int)bins_.size()) std::cout<<"You are trying to access a bin that does not exist!"<<std::endl;
	return &(bins_.at(i-1));
}

unsigned int SearchBins::GetBinNumber(double HT, double MHT, int NJets, int BTags)
{
  unsigned int result =999;
  int match =-1;
  for(unsigned int i=0; i<bins_.size();i++)
  {
	  if(HT>=bins_[i].getHTmin() && 
	      HT<bins_[i].getHTmax() &&
	      MHT>=bins_[i].getMHTmin() && 
	      MHT<bins_[i].getMHTmax() &&
	      NJets+0.1> bins_[i].getNJetsmin() && 
	      NJets-0.1< bins_[i].getNJetsmax() &&
	      BTags+0.1> bins_[i].getBTagsmin() && 
	      BTags-0.1< bins_[i].getBTagsmax()){
		      	result=i;
		      	match++;
		      	usedBin_[i]=usedBin_[i]+1;
	  }
  }
  if(match==-1)
  {
    //std::cout<<"Error event fits in no bin! HT: "<<HT<<", MHT: "<<MHT<<", NJets: "<<NJets<<", BTags: "<<BTags<<std::endl;
  	nSkipped++;
  }
  if(match>0)
  {
    std::cout<<"Error event fits in more than one bin! HT: "<<HT<<", MHT: "<<MHT<<", NJets: "<<NJets<<", BTags: "<<BTags<<std::endl;
  }
  return result+1; // to not use the 0 bin but start counting at 1
}

unsigned int SearchBins::GetCombinedBinNumber(double HT, double MHT, int NJets)
{
  unsigned int result =999;
  int match =-1;
  for(unsigned int i=0; i<binsCombined_.size();i++)
  {
	  if(HT>=binsCombined_[i].getHTmin() && 
	      HT<binsCombined_[i].getHTmax() &&
	      MHT>=binsCombined_[i].getMHTmin() && 
	      MHT<binsCombined_[i].getMHTmax() &&
	      NJets+0.1> binsCombined_[i].getNJetsmin() && 
	      NJets-0.1< binsCombined_[i].getNJetsmax()){
		      	result=i;
		      	match++;
		      	usedCombinedBin_[i]=usedCombinedBin_[i]+1;
	  }
  }
  if(match==-1)
  {
    //std::cout<<"Error event fits in no bin! HT: "<<HT<<", MHT: "<<MHT<<", NJets: "<<NJets<<", BTags: "<<BTags<<std::endl;
  	nSkipped++;
  }
  if(match>0)
  {
    std::cout<<"Error event fits in more than one bin! HT: "<<HT<<", MHT: "<<MHT<<", NJets: "<<NJets<<std::endl;
  }
  return result+1; // to not use the 0 bin but start counting at 1
}

unsigned int SearchBins::GetNbins(){
	return bins_.size();
}

unsigned int SearchBins::GetNbinsCombined(){
	return binsCombined_.size();
}

void SearchBins::PrintUsed()
{
  Long64_t total = 0;
  for(unsigned int i=0; i< usedBin_.size();i++){
  	std::cout<<"Bin["<<i+1<<"] has been used: "<<usedBin_[i]<<std::endl;
  	total+=usedBin_[i];;
  }
  std::cout<<"Total: "<<total<<std::endl;
  std::cout<<"Skipped events (passed baseline but no search bin)/total: "<<nSkipped<<"/"<<nSkipped+total<<std::endl;
}

void SearchBins::PrintUsedCombined()
{
  Long64_t total = 0;
  for(unsigned int i=0; i< usedCombinedBin_.size();i++){
  	std::cout<<"Bin["<<i+1<<"] has been used: "<<usedCombinedBin_[i]<<std::endl;
  	total+=usedCombinedBin_[i];;
  }
  std::cout<<"Total: "<<total<<std::endl;
  std::cout<<"Skipped events (passed baseline but no search bin)/total: "<<nSkipped<<"/"<<nSkipped+total<<std::endl;
}

#endif