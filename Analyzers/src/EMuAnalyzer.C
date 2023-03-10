#include "EMuAnalyzer.h"

EMuAnalyzer::EMuAnalyzer(){
}
EMuAnalyzer::~EMuAnalyzer(){
}
void EMuAnalyzer::executeEvent(){
  //////// Gen level ////////
  if(GetSkimName()==""){
    if(MCSample.Contains("TTLL")){
      FillHist("normcheck",0,reductionweight*MCweight()*_event.GetTriggerLumi("Full"),5,0,5);
      FillHist("normcheck",1,reductionweight*MCweight()*_event.GetTriggerLumi("Full")*mcCorr->GetTopPtReweight(gens),5,0,5);
      FillHist("normcheck",2,reductionweight*MCweight()*_event.GetTriggerLumi("Full")*GetTopPtReweight2(gens),5,0,5);
    }
    return;
  }
  
  //////// nominal channels //////////
  if(!IsDATA||DataStream.Contains("SingleMuon")){
    executeEventWithParameter(MakeParameter("me"));
    for(TString syst:{"jet_scale_up","jet_scale_down","jet_smear_up","jet_smear_down"}){
      executeEventWithParameter(MakeParameter("me",syst));
    }
    //executeEventWithParameter(MakeParameter("me","DeepJet::Medium"));
    //executeEventWithParameter(MakeParameter("me","DeepJet::Tight::mujets"));
    //executeEventWithParameter(MakeParameter("me","DeepCSV::Medium"));
    //executeEventWithParameter(MakeParameter("me","DeepCSV::Tight"));
  }
  if(!IsDATA||DataStream.Contains("SingleElectron")||DataStream.Contains("EGamma")){
    executeEventWithParameter(MakeParameter("em"));
    for(TString syst:{"jet_scale_up","jet_scale_down","jet_smear_up","jet_smear_down"}){
      executeEventWithParameter(MakeParameter("em",syst));
    }
    //executeEventWithParameter(MakeParameter("em","DeepJet::Medium"));
    //executeEventWithParameter(MakeParameter("em","DeepJet::Tight::mujets"));
    //executeEventWithParameter(MakeParameter("em","DeepCSV::Medium"));
    //executeEventWithParameter(MakeParameter("em","DeepCSV::Tight"));
  }
}
SMPAnalyzerCore::Parameter EMuAnalyzer::MakeParameter(TString key,TString option){
  Parameter p=SMPAnalyzerCore::MakeParameter(key,option);
  if(p.suffix==""){
    p.weightbit|=NominalWeight|SystematicWeight|EfficiencyWeight|PDFWeight;
    if(IsTTSample){
      p.weightbit|=PDFWeight;
    }
  }
  if(option.Contains("DeepJet::Medium")){
    p.weightbit=NominalWeight;
    p.prefix="medium/"+p.prefix;
  }else if(option.Contains("DeepJet::Tight::mujets")){
    p.weightbit=NominalWeight;
    p.prefix="mujets/"+p.prefix;
  }else if(option.Contains("DeepCSV::Medium")){
    p.weightbit=NominalWeight;
    p.prefix="csvmedium/"+p.prefix;
  }else if(option.Contains("DeepCSV")){
    p.weightbit=NominalWeight;
    p.prefix="csv/"+p.prefix;
  }
  return p;
}
void EMuAnalyzer::EvalWeights(Parameter& p){
  p.weightmap[""]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.w.btagSF*p.w.topptweight;

  if(!IsDATA){
    p.weightmap["_PUweight_up"]=p.w.lumiweight*p.w.PUweight_up*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.w.btagSF*p.w.topptweight;
    p.weightmap["_PUweight_down"]=p.w.lumiweight*p.w.PUweight_down*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.w.btagSF*p.w.topptweight;

    p.weightmap["_prefireweight_up"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight_up*p.w.zptweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.w.btagSF*p.w.topptweight;
    p.weightmap["_prefireweight_down"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight_down*p.w.zptweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.w.btagSF*p.w.topptweight;

    p.weightmap["_nozptweight"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.w.btagSF*p.w.topptweight;
    p.weightmap["_noz0weight"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.w.btagSF*p.w.topptweight;
    p.weightmap["_noweakweight"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.w.btagSF*p.w.topptweight;
    p.weightmap["_notopptweight"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.w.btagSF;
    p.weightmap["_topptweight2"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.w.btagSF*GetTopPtReweight2(gens);

    p.weightmap["_btagSF_hup"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.w.btagSF_hup;
    p.weightmap["_btagSF_hdown"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.w.btagSF_hdown;
    p.weightmap["_btagSF_lup"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.w.btagSF_lup;
    p.weightmap["_btagSF_ldown"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.w.btagSF_ldown;

    for(int j=0,nj=fEff->nreplica;j<nj;j++){
      double electronRECOSF=p.w.electronRECOSF_sys.size() ? p.w.electronRECOSF_sys[0][j] : 1.;
      double electronIDSF=p.w.electronIDSF_sys.size() ? p.w.electronIDSF_sys[0][j] : 1.;
      double muonIDSF=p.w.muonIDSF_sys.size() ? p.w.muonIDSF_sys[0][j] : 1.;
      double triggerSF=p.w.triggerSF_sys.size() ? p.w.triggerSF_sys[0][j] : 1.;
      p.weightmap[Form("_efficiencySF_stat%d",j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*electronRECOSF*electronIDSF*muonIDSF*p.w.muonISOSF*triggerSF*p.w.CFSF*p.w.btagSF*p.w.topptweight;
    }  
    
    for(int i=1,ni=p.w.electronRECOSF_sys.size();i<ni;i++){
      for(int j=0,nj=p.w.electronRECOSF_sys[i].size();j<nj;j++){
	p.weightmap[Form("_electronRECOSF_s%d_m%d",i,j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF_sys[i][j]*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.w.btagSF*p.w.topptweight;
      }
    }	
    
    for(int i=1,ni=p.w.electronIDSF_sys.size();i<ni;i++){
      for(int j=0,nj=p.w.electronIDSF_sys[i].size();j<nj;j++){
	p.weightmap[Form("_electronIDSF_s%d_m%d",i,j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF_sys[i][j]*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.w.btagSF*p.w.topptweight;
      }
    }	
    
    for(int i=1,ni=p.w.muonIDSF_sys.size();i<ni;i++){
      for(int j=0,nj=p.w.muonIDSF_sys[i].size();j<nj;j++){
	p.weightmap[Form("_muonIDSF_s%d_m%d",i,j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF_sys[i][j]*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.w.btagSF*p.w.topptweight;
      }
    }	
    
    for(int i=1,ni=p.w.triggerSF_sys.size();i<ni;i++){
      for(int j=0,nj=p.w.triggerSF_sys[i].size();j<nj;j++){
	p.weightmap[Form("_triggerSF_s%d_m%d",i,j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF_sys[i][j]*p.w.CFSF*p.w.btagSF*p.w.topptweight;
      }
    }

    if(p.weightbit&PDFWeight){
      for(unsigned int i=0;i<weight_Scale->size();i++){
	p.weightmap[Form("_scalevariation%d",i)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.w.btagSF*p.w.topptweight*weight_Scale->at(i);
      }
      for(unsigned int i=0;i<weight_PDF->size();i++){
	p.weightmap[Form("_pdf%d",i)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.w.btagSF*p.w.topptweight*weight_PDF->at(i);
      }
      if(weight_AlphaS->size()==2){
	p.weightmap["_alphaS_down"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.w.btagSF*p.w.topptweight*weight_AlphaS->at(0);
	p.weightmap["_alphaS_up"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.w.btagSF*p.w.topptweight*weight_AlphaS->at(1);
      }
    }
  }
}
void EMuAnalyzer::FillHists(Parameter& p){
  TLorentzVector dilepton=(*p.lepton0)+(*p.lepton1);
  double dimass=dilepton.M();
  if(dimass<52) return;
  for(const auto& [wname,w]:p.weightmap){
    TString region=p.bjets.size()>0 ? "nbjet/" : "0bjet/";
    TString pre=p.prefix+region+p.hprefix;
    TString suf=p.suffix+wname;
    
    //for leptons
    for(int i=0;i<(int)p.leptons.size();i++){
      double pt=p.leptons.at(i)->Pt();
      double eta=p.leptons.at(i)->Eta();
      FillHist(Form("%sl%dpt%s",pre.Data(),i,suf.Data()),pt,w,nptbin,ptbins);
      FillHist(Form("%sl%deta%s",pre.Data(),i,suf.Data()),eta,w,50,-2.5,2.5);
      
      FillHist(Form("%slpt%s",pre.Data(),suf.Data()),pt,w,nptbin,ptbins);
      FillHist(Form("%sleta%s",pre.Data(),suf.Data()),eta,w,50,-2.5,2.5);

      if(p.leptons.at(i)->LeptonFlavour()==Lepton::Flavour::MUON){
      }else if(p.leptons.at(i)->LeptonFlavour()==Lepton::Flavour::ELECTRON){
	Electron* el=(Electron*)p.leptons.at(i);
	FillHist(Form("%slsceta%s",pre.Data(),suf.Data()),el->scEta(),w,50,-2.5,2.5);
      }
    }

    double dimass=dilepton.M();
    double dipt=dilepton.Pt();
    double dirap=dilepton.Rapidity();
    FillHist(pre+"dimass"+suf,dimass,w,nmbin,mbins);
    FillHist(pre+"dipt"+suf,dipt,w,nptbin,ptbins);
    FillHist(pre+"dirap"+suf,dirap,w,50,-2.5,2.5);
  }
}
