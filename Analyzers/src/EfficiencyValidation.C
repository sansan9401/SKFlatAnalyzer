#include "EfficiencyValidation.h"

EfficiencyValidation::EfficiencyValidation(){
}
EfficiencyValidation::~EfficiencyValidation(){
}
void EfficiencyValidation::initializeAnalyzer(){
  SMPAnalyzerCore::initializeAnalyzer(); //setup zpt roc z0
  fChain->SetBranchStatus("pfMET_*",false);
  fChain->SetBranchStatus("pfMET_Type1_pt",true);
  fChain->SetBranchStatus("jet_*",false);
  fChain->SetBranchStatus("fatjet_*",false);
  fChain->SetBranchStatus("photon_*",false);
}
void EfficiencyValidation::executeEvent(){
  //////// nominal channels //////////
  if(!IsDATA||DataStream.Contains("DoubleMuon")) 
    executeEventWithParameter(MakeParameter("mm"));
  if(!IsDATA||DataStream.Contains("SingleMuon")) 
    executeEventWithParameter(MakeParameter("mu"));
  if(!IsDATA||DataStream.Contains("DoubleEG")||DataStream.Contains("EGamma")) 
    executeEventWithParameter(MakeParameter("ee"));
  if(!IsDATA||DataStream.Contains("SingleElectron")||DataStream.Contains("EGamma")){
    executeEventWithParameter(MakeParameter("el"));
    executeEventWithParameter(MakeParameter("el","TightIDSelectiveCharge"));
  }

  //////// testing channels //////////
  
  if(GetEra()=="2017"){
    if(!IsDATA||DataStream.Contains("DoubleEG")){
      Parameter p=MakeParameter("ee");
      p.suffix="_noL1";
      p.k.triggerSF={"Ele23Leg1_MediumID_Q_v3_2","Ele12Leg2_MediumID_Q"};
      //executeEventWithParameter(p);
    }
  }else if(GetEra()=="2018"){
    if(!IsDATA||DataStream.Contains("DoubleMuon")){
      Parameter p=MakeParameter("mm");
      p.suffix="_new";
      p.SetMuonKeys("IDISO_SF_MediumID_trkIsoLoose_Q_new","",{"Mu17Leg1_MediumID_trkIsoLoose_Q_new","Mu8Leg2_MediumID_trkIsoLoose_Q_new"});
      //executeEventWithParameter(p);
    }
    if(!IsDATA||DataStream.Contains("SingleMuon")){ 
      Parameter p=MakeParameter("mu");
      p.suffix="_new";
      p.SetMuonKeys("IDISO_SF_MediumID_trkIsoLoose_Q_new","",{"IsoMu24_MediumID_trkIsoLoose_Q_new"});
      //executeEventWithParameter(p);
    }
    if(!IsDATA||DataStream.Contains("EGamma")){
      Parameter p=MakeParameter("ee");
      p.suffix="_noL1";
      p.k.triggerSF={"Ele23Leg1_MediumID_Q_v3","Ele12Leg2_MediumID_Q"};
      //executeEventWithParameter(p);

      p=MakeParameter("el");
      p.prefix="el201828/";
      p.triggers={"HLT_Ele28_WPTight_Gsf_v"};
      p.k.triggerSF={"Ele28_MediumID_Q"};
      if(!IsDATA) p.w.lumiweight*=23687.253/59827.879;
      //executeEventWithParameter(p);

      p=MakeParameter("el");
      p.prefix="el201832/";
      p.triggers={"HLT_Ele32_WPTight_Gsf_v"};
      p.k.triggerSF={"Ele32_MediumID_Q"};
      p.c.lepton0pt=35;
      //executeEventWithParameter(p);
    }
  }    
}
SMPAnalyzerCore::Parameter EfficiencyValidation::MakeParameter(TString key,TString option){
  Parameter p=SMPAnalyzerCore::MakeParameter(key,option);
  p.weightbit|=EfficiencyWeight;
  return p;
}
void EfficiencyValidation::FillHists(Parameter& p){
  p.weightmap[""]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
  if(!IsDATA&&p.suffix==""){
    p.weightmap["_noweight"]=p.w.lumiweight;
    p.weightmap["_noPUweight"]=p.w.lumiweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
    p.weightmap["_noprefireweight"]=p.w.lumiweight*p.w.PUweight*p.w.zptweight*p.w.z0weight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
    p.weightmap["_nozptweight"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
    p.weightmap["_noz0weight"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
    
    p.weightmap["_noelectronRECOSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
    for(int i=0;i<p.w.electronRECOSF_sys.size();i++){
      for(int j=0;j<p.w.electronRECOSF_sys[i].size();j++){
	p.weightmap[Form("_electronRECOSF_s%d_m%d",i,j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.electronRECOSF_sys[i][j]*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
      }
    }	
    
    p.weightmap["_noelectronIDSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.electronRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
    for(int i=0;i<p.w.electronIDSF_sys.size();i++){
      for(int j=0;j<p.w.electronIDSF_sys[i].size();j++){
	p.weightmap[Form("_electronIDSF_s%d_m%d",i,j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.electronRECOSF*p.w.electronIDSF_sys[i][j]*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
      }
    }	

    p.weightmap["_nomuonIDSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
    for(int i=0;i<p.w.muonIDSF_sys.size();i++){
      for(int j=0;j<p.w.muonIDSF_sys[i].size();j++){
	p.weightmap[Form("_muonIDSF_s%d_m%d",i,j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF_sys[i][j]*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
      }
    }	
    
    p.weightmap["_noISOSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.electronRECOSF*p.w.electronIDSF*p.w.triggerSF*p.w.CFSF;
    for(int i=0;i<p.w.muonISOSF_sys.size();i++){
      for(int j=0;j<p.w.muonISOSF_sys[i].size();j++){
	p.weightmap[Form("_ISOSF_s%d_m%d",i,j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF_sys[i][j]*p.w.triggerSF*p.w.CFSF;
      }
    }	
    
    p.weightmap["_notriggerSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.CFSF;
    for(int i=0;i<p.w.triggerSF_sys.size();i++){
      for(int j=0;j<p.w.triggerSF_sys[i].size();j++){
	p.weightmap[Form("_triggerSF_s%d_m%d",i,j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF_sys[i][j]*p.w.CFSF;
      }
    }	
    
    p.weightmap["_noefficiencySF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.CFSF;
    
    p.weightmap["_noCFSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF;
    p.weightmap["_CFSF_up"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF_up;
    p.weightmap["_CFSF_down"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF_down;
  }

  TLorentzVector dilepton=(*p.lepton0)+(*p.lepton1);
  double dimass=dilepton.M();
  if(dimass>=52&&dimass<150){
    FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"m52to150",p.weightmap[""]);
    FillHistsEfficiency(p,"m52to150/");
  }
  if(dimass>=80&&dimass<100){
    FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"m80to100",p.weightmap[""]);
    FillHistsEfficiency(p,"m80to100/");
  }
}

void EfficiencyValidation::FillHistsEfficiency(Parameter& p,TString region){
  TString pre=p.prefix+region+p.hprefix;
  for(const auto& [wname,w]:p.weightmap){
    TString suf=p.suffix+wname;
    
    //for leptons
    for(int i=0;i<(int)p.leptons.size();i++){
      double pt=p.leptons.at(i)->Pt();
      double eta=p.leptons.at(i)->Eta();
      TString charge=p.leptons.at(i)->Charge()>0?"p":"m";
      FillHist(Form("%sl%dpt%s",pre.Data(),i,suf.Data()),pt,w,500,0,500);
      FillHist(Form("%sl%deta%s",pre.Data(),i,suf.Data()),eta,w,120,-3,3);
      
      FillHist(Form("%slpt%s",pre.Data(),suf.Data()),pt,w,500,0,500);
      FillHist(Form("%sleta%s",pre.Data(),suf.Data()),eta,w,120,-3,3);
    }

    TLorentzVector dilepton=(*p.lepton0)+(*p.lepton1);
    double dimass=dilepton.M();
    double dipt=dilepton.Pt();
    double dirap=dilepton.Rapidity();
    FillHist(pre+"dimass"+suf,dimass,w,196,52,150);
    FillHist(pre+"dipt"+suf,dipt,w,400,0,400);
    FillHist(Form("%snPV%s",pre.Data(),suf.Data()),nPV,w,100,0,100);
    if(!IsDATA) FillHist(Form("%snPileUp%s",pre.Data(),suf.Data()),nPileUp,w,100,0,100);


    if(wname!="") continue;

    //for leptons
    for(int i=0;i<(int)p.leptons.size();i++){
      if(i>1) break;
      double pt=p.leptons.at(i)->Pt();
      double eta=p.leptons.at(i)->Eta();
      TString charge=p.leptons.at(i)->Charge()>0?"p":"m";

      FillHist(Form("%sl%d%spt%s",pre.Data(),i,charge.Data(),suf.Data()),pt,w,500,0,500);
      FillHist(Form("%sl%d%seta%s",pre.Data(),i,charge.Data(),suf.Data()),eta,w,120,-3,3);
      
      FillHist(Form("%sl%spt%s",pre.Data(),charge.Data(),suf.Data()),pt,w,500,0,500);
      FillHist(Form("%sl%seta%s",pre.Data(),charge.Data(),suf.Data()),eta,w,120,-3,3);
      
      FillHist(Form("%slriso%s",pre.Data(),suf.Data()),p.leptons.at(i)->RelIso(),w,30,0,0.3);
      if(p.leptons.at(i)->LeptonFlavour()==Lepton::Flavour::MUON){
	double rtrkiso=((Muon*)p.leptons.at(i))->TrkIso()/pt;
	FillHist(Form("%slrtrkiso%s",pre.Data(),suf.Data()),rtrkiso,w,40,0,0.2);
      }else if(p.leptons.at(i)->LeptonFlavour()==Lepton::Flavour::ELECTRON){
	Electron* el=(Electron*)p.leptons.at(i);
	FillHist(Form("%slrawpt%s",pre.Data(),suf.Data()),el->UncorrPt(),w,500,0,500);
	FillHist(Form("%slsceta%s",pre.Data(),suf.Data()),el->scEta(),w,120,-3,3);
      }
    }
      
    FillHist(pre+"dirap"+suf,dirap,w,120,-3,3);
    FillHist(pre+"nlepton"+suf,p.muons.size()+p.electrons.size(),w,10,0,10);
    FillHist(pre+"met"+suf,pfMET_Type1_pt,w,100,0,200);
  }
}
