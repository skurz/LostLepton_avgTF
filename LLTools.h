#ifndef LLTOOLS_H
#define LLTOOLS_H

#include "TVector2.h"
#include "TH2.h"
#include <cmath>
#include <iostream>
#include <utility>
#include <TROOT.h>
#include <TCanvas.h>

static double deltaR(double eta1, double phi1, double eta2, double phi2)
{
  double deta = eta1-eta2;
  double dphi = TVector2::Phi_mpi_pi(phi1-phi2);
  return sqrt(deta * deta + dphi *dphi); 
}

static double MTWCalculator(double metPt,double  metPhi,double  lepPt,double  lepPhi, int scaleMet = 0)
{
  if(scaleMet == 1) metPt = metPt * 1.3;
  if(scaleMet == -1) metPt = metPt * 0.7;
  double deltaPhi =TVector2::Phi_mpi_pi(lepPhi-metPhi);
  return sqrt(2*lepPt*metPt*(1-cos(deltaPhi)) );
}

static double PTWCalculator(double metPt,double  metPhi,double  lepPt,double  lepPhi)
{
  TVector2 Wvec(metPt*cos(metPhi)+lepPt*cos(lepPhi), metPt*sin(metPhi)+lepPt*sin(lepPhi));
  return Wvec.Mod();
}

static double GetCosDTT(double mht, double mht_phi, double lep_pt, double lep_phi) {

  double px_lab_w = lep_pt * cos( lep_phi ) + mht*cos(mht_phi) ;
  double py_lab_w = lep_pt * sin( lep_phi ) + mht*sin(mht_phi) ;
  double pt_lab_w = sqrt( pow( px_lab_w, 2 ) + pow( py_lab_w, 2 ) ) ;

  double dphi_lab_lepton_mht = TVector2::Phi_mpi_pi(lep_phi - mht_phi);
  double cosdphi_lab_lepton_mht = cos( dphi_lab_lepton_mht ) ;

  double mt_w = sqrt( 2*lep_pt*mht*(1.-cosdphi_lab_lepton_mht) ) ;
  double et_lab_w = sqrt( pow( pt_lab_w, 2. ) + pow( mt_w, 2. ) ) ;

	 
  double pt_parallel_lab_lepton = (lep_pt/pt_lab_w)*(lep_pt + mht*cosdphi_lab_lepton_mht ) ;

  return (2*pt_parallel_lab_lepton - pt_lab_w)/et_lab_w ;
  
}

static double GetTriggerEffWeight(Double_t MHT){
  if (MHT<50) return 0.0080;
  else if (MHT<100) return 0.0988;
  else if (MHT<150) return 0.5424;
  else if (MHT<200) return 0.9011;
  else if (MHT<250) return 0.9846;
  else if (MHT<300) return 0.9961;
  else if (MHT<400) return 0.9865;
  else if (MHT<500) return 1.0000;
  else return 1.0000;
}

static double GetSignalTriggerEffWeight(Double_t MHT) {
  if (MHT<250) return 0;
  else if(MHT<300) return 0.982;
  else if (MHT<350) return 0.985;
  else if (MHT<500) return 0.995;
  else return 1.00;
}

static std::pair<double,double> EvalSF(TH2 *hist, Double_t xVal, Double_t yVal) {
  // Dont use overflow bins!
  if(xVal < hist->GetXaxis()->GetXmin() )
    {
      //std::cout<<"SF: Warning xVal: "<<xVal<<" is smaller than minimum of histo: "<<hist->GetName()<<std::endl;
      xVal= hist->GetXaxis()->GetXmin()+0.01;
    }
  else if(xVal > hist->GetXaxis()->GetXmax() )
    {
      //std::cout<<"SF: Warning xVal: "<<xVal<<" is bigger than maximum of histo: "<<hist->GetName()<<" which is: "<<hist->GetXaxis()->GetXmax()<<std::endl;
      xVal= hist->GetXaxis()->GetXmax()-0.01;
    }
  
  if(yVal < hist->GetYaxis()->GetXmin() )
    {
      //std::cout<<"SF: Warning yVal: "<<yVal<<" is smaller than minimum of histo: "<<hist->GetName()<<std::endl;
      yVal= hist->GetYaxis()->GetXmin()+0.01;
    }
  else if(yVal > hist->GetYaxis()->GetXmax() )
    {
      //std::cout<<"SF: Warning yVal: "<<yVal<<" is bigger than maximum of histo: "<<hist->GetName()<<std::endl;
      yVal= hist->GetYaxis()->GetXmax()-0.01;
    }

  int nxBin = hist->GetXaxis()->FindBin(xVal);
  int nyBin = hist->GetYaxis()->FindBin(yVal);

  if(nxBin > hist->GetNbinsX() || nyBin > hist->GetNbinsY()) std::cout<<"SF: Problem in getting Efficiencies!"<<std::endl;
  if(nxBin > hist->GetNbinsX()) nxBin = hist->GetNbinsX();
  if(nyBin > hist->GetNbinsY()) nyBin = hist->GetNbinsY();

  return std::make_pair(hist->GetBinContent(nxBin, nyBin), hist->GetBinError(nxBin, nyBin));
}

static double GetSFUnc(TH2 *hist, Double_t xVal, Double_t yVal, double addSys = 0.) {

  std::pair<double, double> SFandUnc = EvalSF(hist, xVal, yVal);

  double SF = 0.;

  if(addSys > 0) SF = std::sqrt(SFandUnc.second*SFandUnc.second + addSys*addSys);
  else SF = SFandUnc.second;

  //std::cout << std::abs(1-hist->GetBinContent(nxBin, nyBin)) << " " << std::sqrt(hist->GetBinError(nxBin, nyBin)*hist->GetBinError(nxBin, nyBin) + 0.01*hist->GetBinContent(nxBin, nyBin)*0.01*hist->GetBinContent(nxBin, nyBin)) << " " << hist->GetBinError(nxBin, nyBin)<<std::endl;

  return SF;
}

static double GetSF(TH2 *hist, Double_t xVal, Double_t yVal) {
  return EvalSF(hist, xVal, yVal).first;
}

static std::pair<double,double> EvalSF(TH1 *hist, Double_t xVal) {
  // Dont use overflow bins!
  if(xVal < hist->GetXaxis()->GetXmin() )
    {
      //std::cout<<"SF: Warning xVal: "<<xVal<<" is smaller than minimum of histo: "<<hist->GetName()<<std::endl;
      xVal= hist->GetXaxis()->GetXmin()+0.01;
    }
  else if(xVal > hist->GetXaxis()->GetXmax() )
    {
      //std::cout<<"SF: Warning xVal: "<<xVal<<" is bigger than maximum of histo: "<<hist->GetName()<<" which is: "<<hist->GetXaxis()->GetXmax()<<std::endl;
      xVal= hist->GetXaxis()->GetXmax()-0.01;
    }
  
  int nxBin = hist->GetXaxis()->FindBin(xVal);

  if(nxBin > hist->GetNbinsX()) std::cout<<"SF: Problem in getting Efficiencies!"<<std::endl;
  if(nxBin > hist->GetNbinsX()) nxBin = hist->GetNbinsX();

  return std::make_pair(hist->GetBinContent(nxBin), hist->GetBinError(nxBin));
}

static double GetSFUnc(TH1 *hist, Double_t xVal, double addSys = 0.) {

  std::pair<double, double> SFandUnc = EvalSF(hist, xVal);

  double SF = 0.;

  if(addSys > 0) SF = std::sqrt(SFandUnc.second*SFandUnc.second + addSys*addSys);
  else SF = SFandUnc.second;

  return SF;
}

static double GetSF(TH1 *hist, Double_t xVal) {
  return EvalSF(hist, xVal).first;
}

static double GetAccSF(std::vector<TH2D*> hists, Int_t NJets, Double_t xVal, Double_t yVal) {
  if(hists.size() != 7){
    std::cout<<"Wrong dimension of histograms for acceptance uncertainty maps!"<<std::endl;
    return 0;
  }

  TH2D* hist = 0;

  switch(NJets){
    case 2:
      hist = hists.at(0);
      break;
    case 3:
      hist = hists.at(1);
      break;
    case 4:
      hist = hists.at(2);
      break;
    case 5:
      hist = hists.at(3);
      break;
    case 6:
      hist = hists.at(4);
      break;
    case 7:
    case 8:
      hist = hists.at(5);
      break;
    default:
      hist = hists.at(6);
      break;
  }

  // Dont use overflow bins!
  if(xVal < hist->GetXaxis()->GetXmin() )
    {
      //std::cout<<"SF: Warning xVal: "<<xVal<<" is smaller than minimum of histo: "<<hist->GetName()<<std::endl;
      xVal= hist->GetXaxis()->GetXmin()+0.01;
    }
  else if(xVal > hist->GetXaxis()->GetXmax() )
    {
      //std::cout<<"SF: Warning xVal: "<<xVal<<" is bigger than maximum of histo: "<<hist->GetName()<<" which is: "<<hist->GetXaxis()->GetXmax()<<std::endl;
      xVal= hist->GetXaxis()->GetXmax()-0.01;
    }
  
  if(yVal < hist->GetYaxis()->GetXmin() )
    {
      //std::cout<<"SF: Warning yVal: "<<yVal<<" is smaller than minimum of histo: "<<hist->GetName()<<std::endl;
      yVal= hist->GetYaxis()->GetXmin()+0.01;
    }
  else if(yVal > hist->GetYaxis()->GetXmax() )
    {
      //std::cout<<"SF: Warning yVal: "<<yVal<<" is bigger than maximum of histo: "<<hist->GetName()<<std::endl;
      yVal= hist->GetYaxis()->GetXmax()-0.01;
    }

  int nxBin = hist->GetXaxis()->FindBin(xVal);
  int nyBin = hist->GetYaxis()->FindBin(yVal);

  if(nxBin > hist->GetNbinsX() || nyBin > hist->GetNbinsY()) std::cout<<"SF: Problem in getting Efficiencies!"<<std::endl;
  if(nxBin > hist->GetNbinsX()) nxBin = hist->GetNbinsX();
  if(nyBin > hist->GetNbinsY()) nyBin = hist->GetNbinsY();

  return hist->GetBinContent(nxBin, nyBin);
}

static void savePlot(TH1* plot){
  gROOT->SetBatch(true);
  TString name_(plot->GetName());
  TCanvas *c1 = new TCanvas(name_,plot->GetTitle(),1);
  c1->cd();
/*  if (xlog) {
    c1->SetLogx();
    RatioTH1D_->GetXaxis()->SetRangeUser(0.001, 10*RatioTH1D_->GetBinLowEdge(RatioTH1D_->GetNbinsX()+1));
  }
  if (ylog) c1->SetLogy();
  if(RatioTH1D_->GetXaxis()->GetBinCenter(RatioTH1D_->GetXaxis()->GetNbins()) > 500) c1->SetLogx();
  RatioTH1D_->GetYaxis()->SetRangeUser(0.001, 1.05);
*/
  plot->Draw("Text");
  c1->SaveAs("AdditionalPlots/"+name_+".pdf");
  delete c1;
  gROOT->SetBatch(false);
}


static double getMuonIDSF(Double_t pt, Double_t eta){
  double sf = 1.;

  if(std::abs(eta)<1.2){
    if(pt < 20) sf = 0.983;
    else if(pt < 30) sf = 0.985;
    else if(pt < 40) sf = 0.992;
    else if(pt < 50) sf = 0.994;
    else if(pt < 60) sf = 0.988;
    else if(pt < 80) sf = 0.986;
    else if(pt < 120) sf = 0.979;
    else sf = 1.016;
  }
  else{
    if(pt < 20) sf = 1.004;
    else if(pt < 30) sf = 0.994;
    else if(pt < 40) sf = 0.99;
    else if(pt < 50) sf = 0.992;
    else if(pt < 60) sf = 0.986;
    else if(pt < 80) sf = 0.981;
    else if(pt < 120) sf = 0.955;
    else sf = 1.005;
  }

  return std::abs(1. - sf);
}

static double getMuonIsoSF(Double_t pt, Double_t eta, Double_t act){
  double sf = 1.;

  if(std::abs(eta)<1.2){
    if(pt < 20) sf = 1.002;
    else if(pt < 30) sf = 1.001;
    else if(pt < 40) sf = 1.000;
    else if(pt < 50) sf = 1.000;
    else if(pt < 60) sf = 1.001;
    else if(pt < 80) sf = 1.001;
    else if(pt < 120) sf = 1.000;
    else sf = 1.000;
  }
  else{
    if(pt < 20) sf = 0.995;
    else if(pt < 30) sf = 1.001;
    else if(pt < 40) sf = 1.000;
    else if(pt < 50) sf = 1.000;
    else if(pt < 60) sf = 1.000;
    else if(pt < 80) sf = 1.000;
    else if(pt < 120) sf = 1.000;
    else sf = 1.000;
  }

  return std::abs(1. - sf);
}

static double getElecIDSF(Double_t pt, Double_t eta){
  double sf = 1.;

  if(std::abs(eta)<1.442){
    if(pt < 20) sf = 0.941;
    else if(pt < 30) sf = 0.981;
    else if(pt < 40) sf = 0.995;
    else if(pt < 50) sf = 0.996;
    else sf = 1.001;
  }else if(std::abs(eta)<1.5){
    if(pt < 20) sf = 1.135;
    else if(pt < 30) sf = 1.035;
    else if(pt < 40) sf = 0.974;
    else if(pt < 50) sf = 1.006;
    else sf = 1.000;
  }else{
    if(pt < 20) sf = 0.970;
    else if(pt < 30) sf = 0.965;
    else if(pt < 40) sf = 1.002;
    else if(pt < 50) sf = 1.016;
    else sf = 1.006;
  }

  return std::abs(1. - sf);
}

static double getElecIsoSF(Double_t pt, Double_t eta, Double_t act){
  double sf = 1.;

  if(act<0.02){
    if(pt < 20) sf = 0.988;
    else if(pt < 30) sf = 0.996;
    else if(pt < 40) sf = 0.996;
    else if(pt < 50) sf = 0.997;
    else sf = 1.001;
  }else if(act<0.05){
    if(pt < 20) sf = 0.948;
    else if(pt < 30) sf = 0.991;
    else if(pt < 40) sf = 0.997;
    else if(pt < 50) sf = 0.996;
    else sf = 0.995;
  }else if(act<0.15){
    if(pt < 20) sf = 1.007;
    else if(pt < 30) sf = 0.992;
    else if(pt < 40) sf = 0.997;
    else if(pt < 50) sf = 0.997;
    else sf = 0.996;
  }else if(act<1.){
    if(pt < 20) sf = 0.992;
    else if(pt < 30) sf = 0.972;
    else if(pt < 40) sf = 0.989;
    else if(pt < 50) sf = 1.004;
    else sf = 0.977;
  }else{
    if(pt < 20) sf = 0.826;
    else if(pt < 30) sf = 0.952;
    else if(pt < 40) sf = 1.090;
    else if(pt < 50) sf = 0.891;
    else sf = 0.930;
  }

  return std::abs(1. - sf);
}

#endif