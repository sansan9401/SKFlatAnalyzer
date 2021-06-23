#include "ZpeakAnalyzer.h"

ZpeakAnalyzer::ZpeakAnalyzer(){
}
ZpeakAnalyzer::~ZpeakAnalyzer(){
}
void ZpeakAnalyzer::initializeAnalyzer(){
  SMPAnalyzerCore::initializeAnalyzer(); //setup zpt roc z0
  fChain->SetBranchStatus("pfMET_*",false);
  fChain->SetBranchStatus("pfMET_Type1_pt",true);
  fChain->SetBranchStatus("jet_*",false);
  fChain->SetBranchStatus("fatjet_*",false);
  fChain->SetBranchStatus("photon_*",false);
}
void ZpeakAnalyzer::executeEvent(){
  if(!IsDATA||DataStream.Contains("DoubleMuon")){
    executeEventWithParameter(MakeParameter("mm"));
    executeEventWithParameter(MakeParameter("MM"));
  }
  if(!IsDATA||DataStream.Contains("DoubleEG")||DataStream.Contains("EGamma")){
    executeEventWithParameter(MakeParameter("ee"));
    executeEventWithParameter(MakeParameter("EE"));
  }
}

void ZpeakAnalyzer::FillHists(Parameter& p){
  map<TString,double> weightmap;
  weightmap[""]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
  weightmap["_noCFSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF;

  for(auto [suf,weight]:weightmap){
    TLorentzVector dilepton=(*p.lepton0)+(*p.lepton1);
    double dimass=dilepton.M();
    FillHist(p.prefix+p.hprefix+"dimass"+p.suffix+suf,p.lepton0->Pt(),p.lepton0->Eta(),dimass,weight,nptbin,ptbins,netabin,etabins,nmassbin,massbins);
    FillHist(p.prefix+p.hprefix+"dimass"+p.suffix+suf,p.lepton1->Pt(),p.lepton1->Eta(),dimass,weight,nptbin,ptbins,netabin,etabins,nmassbin,massbins);
    TH1* h=GetHist3D(p.prefix+p.hprefix+"dimass"+p.suffix+suf);
    int l0ptbin=h->GetXaxis()->FindBin(p.lepton0->Pt());
    int l0etabin=h->GetYaxis()->FindBin(p.lepton0->Eta());
    int l1ptbin=h->GetXaxis()->FindBin(p.lepton1->Pt());
    int l1etabin=h->GetYaxis()->FindBin(p.lepton1->Eta());
    if((l0ptbin-l1ptbin+l0etabin-l1etabin)%2==0){
      FillHist(p.prefix+p.hprefix+"dimass_even"+p.suffix+suf,p.lepton0->Pt(),p.lepton0->Eta(),dimass,weight,nptbin,ptbins,netabin,etabins,nmassbin,massbins);
      FillHist(p.prefix+p.hprefix+"dimass_even"+p.suffix+suf,p.lepton1->Pt(),p.lepton1->Eta(),dimass,weight,nptbin,ptbins,netabin,etabins,nmassbin,massbins);
    }else{
      FillHist(p.prefix+p.hprefix+"dimass_odd"+p.suffix+suf,p.lepton0->Pt(),p.lepton0->Eta(),dimass,weight,nptbin,ptbins,netabin,etabins,nmassbin,massbins);
      FillHist(p.prefix+p.hprefix+"dimass_odd"+p.suffix+suf,p.lepton1->Pt(),p.lepton1->Eta(),dimass,weight,nptbin,ptbins,netabin,etabins,nmassbin,massbins);
    }
    if(!IsDATA){
      if(p.lepton0->Charge()*p.truth_lepton0.Charge()<0){
	FillHist(p.prefix+p.hprefix+"dimass_cf"+p.suffix+suf,p.lepton0->Pt(),p.lepton0->Eta(),dimass,weight,nptbin,ptbins,netabin,etabins,nmassbin,massbins);
	if((l0ptbin-l1ptbin+l0etabin-l1etabin)%2==0)
	  FillHist(p.prefix+p.hprefix+"dimass_cf_even"+p.suffix+suf,p.lepton0->Pt(),p.lepton0->Eta(),dimass,weight,nptbin,ptbins,netabin,etabins,nmassbin,massbins);
	else
	  FillHist(p.prefix+p.hprefix+"dimass_cf_odd"+p.suffix+suf,p.lepton0->Pt(),p.lepton0->Eta(),dimass,weight,nptbin,ptbins,netabin,etabins,nmassbin,massbins);
      }
      if(p.lepton1->Charge()*p.truth_lepton1.Charge()<0){
	FillHist(p.prefix+p.hprefix+"dimass_cf"+p.suffix+suf,p.lepton1->Pt(),p.lepton1->Eta(),dimass,weight,nptbin,ptbins,netabin,etabins,nmassbin,massbins);
	if((l0ptbin-l1ptbin+l0etabin-l1etabin)%2==0)
	  FillHist(p.prefix+p.hprefix+"dimass_cf_even"+p.suffix+suf,p.lepton1->Pt(),p.lepton1->Eta(),dimass,weight,nptbin,ptbins,netabin,etabins,nmassbin,massbins);
	else
	  FillHist(p.prefix+p.hprefix+"dimass_cf_odd"+p.suffix+suf,p.lepton1->Pt(),p.lepton1->Eta(),dimass,weight,nptbin,ptbins,netabin,etabins,nmassbin,massbins);
      }
    }
  }
}
