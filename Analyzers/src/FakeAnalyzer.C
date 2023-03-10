#include "FakeAnalyzer.h"

FakeAnalyzer::FakeAnalyzer(){
}
FakeAnalyzer::~FakeAnalyzer(){
}
void FakeAnalyzer::initializeAnalyzer(){
  SMPAnalyzerCore::initializeAnalyzer(); //setup zpt roc z0
  jtp = JetTagging::Parameters(JetTagging::DeepJet,JetTagging::Tight,JetTagging::incl,JetTagging::comb);
  vector<JetTagging::Parameters> jtps={jtp};
  mcCorr->SetJetTaggingParameters(jtps);
}
void FakeAnalyzer::executeEvent(){
  if(GetSkimName()=="Dilepton") executeDileptonEvent();
  else if(GetSkimName()=="") executeLeptonEvent();
}
void FakeAnalyzer::executeLeptonEvent(){
  if(!IsDATA
     || (DataStream.Contains("DoubleEG")&&DataYear==2016)
     || (DataStream.Contains("SingleElectron")&&DataYear==2017)
     || (DataStream.Contains("EGamma")&&DataYear==2018) ){

    executeEventWithParameter(MakeParameter("ej"));    
    executeEventWithParameter(MakeParameter("Ej"));
  }
  if(!IsDATA||DataStream.Contains("DoubleMuon")){
    executeEventWithParameter(MakeParameter("mj"));
    executeEventWithParameter(MakeParameter("Mj"));
  }
}
void FakeAnalyzer::executeDileptonEvent(){
  if(!IsDATA||DataStream.Contains("DoubleMuon")){
    {Parameter p=MakeParameter("mm");
    p.SetAMuons(MuonMomentumCorrection(SMPGetMuons("POGMediumWithAntiLooseTrkIso",0.0,2.4),0,0));
    p.option+=" triggermatching strictorder";
    executeEventWithParameter(p);}

    {Parameter p=MakeParameter("mm","nobjetcleaning");
    p.SetAMuons(MuonMomentumCorrection(SMPGetMuons("POGMediumWithAntiLooseTrkIso",0.0,2.4),0,0));
    p.option+=" triggermatching strictorder";
    p.prefix="nobjetcleaning/"+p.prefix;
    executeEventWithParameter(p);}
      
    executeEventWithParameter(MakeParameter("mM"));
    executeEventWithParameter(MakeParameter("Mm"));
    executeEventWithParameter(MakeParameter("MM"));

    {Parameter p=MakeParameter("mM");
    p.SetAMuons(ToConePt(p.amuons));
    p.prefix+="cpt/";
    executeEventWithParameter(p);}

    {Parameter p=MakeParameter("Mm");
    p.SetAMuons(ToConePt(p.amuons));
    p.prefix+="cpt/";
    executeEventWithParameter(p);}

    {Parameter p=MakeParameter("MM");
    p.SetAMuons(ToConePt(p.amuons));
    p.prefix+="cpt/";
    executeEventWithParameter(p);}

    {Parameter p=MakeParameter("mM");
    p.SetAMuons(ToModifiedPt(p.amuons));
    p.prefix+="mpt/";
    executeEventWithParameter(p);}

    {Parameter p=MakeParameter("Mm");
    p.SetAMuons(ToModifiedPt(p.amuons));
    p.prefix+="mpt/";
    executeEventWithParameter(p);}

    {Parameter p=MakeParameter("MM");
    p.SetAMuons(ToModifiedPt(p.amuons));
    p.prefix+="mpt/";
    executeEventWithParameter(p);}

  }
  if(!IsDATA||DataStream.Contains("DoubleEG")||DataStream.Contains("EGamma")){
    {Parameter p=MakeParameter("ee");
    p.SetAElectrons(SMPGetElectrons("passMediumIDSideBand",0.0,2.5));
    p.option+=" triggermatching strictorder";
    executeEventWithParameter(p);}

    {Parameter p=MakeParameter("ee","nobjetcleaning");
    p.SetAElectrons(SMPGetElectrons("passMediumIDSideBand",0.0,2.5));
    p.option+=" triggermatching strictorder";
    p.prefix="nobjetcleaning/"+p.prefix;
    executeEventWithParameter(p);}

    executeEventWithParameter(MakeParameter("eE"));
    executeEventWithParameter(MakeParameter("Ee"));
    executeEventWithParameter(MakeParameter("EE"));

    {Parameter p=MakeParameter("eE");
    p.SetAElectrons(ToConePt(p.aelectrons));
    p.prefix+="cpt/";
    executeEventWithParameter(p);}

    {Parameter p=MakeParameter("Ee");
    p.SetAElectrons(ToConePt(p.aelectrons));
    p.prefix+="cpt/";
    executeEventWithParameter(p);}

    {Parameter p=MakeParameter("EE");
    p.SetAElectrons(ToConePt(p.aelectrons));
    p.prefix+="cpt/";
    executeEventWithParameter(p);}

    {Parameter p=MakeParameter("eE");
    p.SetAElectrons(ToModifiedPt(p.aelectrons));
    p.prefix+="mpt/";
    executeEventWithParameter(p);}

    {Parameter p=MakeParameter("Ee");
    p.SetAElectrons(ToModifiedPt(p.aelectrons));
    p.prefix+="mpt/";
    executeEventWithParameter(p);}

    {Parameter p=MakeParameter("EE");
    p.SetAElectrons(ToModifiedPt(p.aelectrons));
    p.prefix+="mpt/";
    executeEventWithParameter(p);}
  }
}

void FakeAnalyzer::FillHists(Parameter& p){
  if(GetSkimName()=="Dilepton"){
    TLorentzVector dilepton=(*p.lepton0)+(*p.lepton1);
    double dimass=dilepton.M();
    if(dimass>52){
      if(p.bjets.size()) FillDileptonHists(p,"nbjet/");
      else FillDileptonHists(p,"0bjet/");
      FillDileptonHists(p,"");
      if(dimass<76||dimass>106){
	if(p.bjets.size()) FillDileptonHists(p,"nbjet/noZ/");
	else FillDileptonHists(p,"0bjet/noZ/");
	FillDileptonHists(p,"noZ/");
      }
    } 
  }else if(GetSkimName()==""){
    std::vector<Jet> jets30=GetJets("tightLepVeto",30,2.4);
    std::sort(jets30.begin(),jets30.end(),PtComparing);
    if(jets30.size()&&fabs(jets30.at(0).DeltaPhi(*p.lepton0))>TMath::Pi()*2/3&&(jets30.at(0)+*p.lepton0).M()>50){
      if(p.bjets.size()) FillLeptonHists(p,"nbjet/");
      else FillLeptonHists(p,"0bjet/");
      FillLeptonHists(p,"");
    }
  }
}
void FakeAnalyzer::FillLeptonHists(Parameter& p,TString region){
  for(auto [suf,weight]:p.weightmap){
    double eta=fabs(p.leptons.at(0)->Eta());
    double pt=p.leptons.at(0)->Pt();
    double riso=p.leptons.at(0)->RelIso();
    double energy=p.leptons.at(0)->E();
    FillHist(p.prefix+region+p.hprefix+"l0pt"+p.suffix+suf,pt,weight,nptbin,ptbins);
    FillHist(p.prefix+region+p.hprefix+"l0eta"+p.suffix+suf,eta,weight,netabin,etabins);
    FillHist(p.prefix+region+p.hprefix+"l0etapt"+p.suffix+suf,eta,pt,weight,netabin,etabins,nptbin,ptbins);
    FillHist(p.prefix+region+p.hprefix+"l0etae"+p.suffix+suf,eta,energy,weight,netabin,etabins,nptbin,ptbins);
    FillHist(p.prefix+region+p.hprefix+"l0riso"+p.suffix+suf,riso,weight,100,0,2);
  }
}
void FakeAnalyzer::FillDileptonHists(Parameter& p,TString region){
  for(auto [suf,weight]:p.weightmap){
    TLorentzVector dilepton=(*p.lepton0)+(*p.lepton1);
    FillHist(p.prefix+region+p.hprefix+"dimass"+p.suffix+suf,dilepton.M(),weight,98,52,150);
    FillHist(p.prefix+region+p.hprefix+"dimass_wide"+p.suffix+suf,dilepton.M(),weight,nmbin,mbins);
    FillHist(p.prefix+region+p.hprefix+"dipt"+p.suffix+suf,dilepton.Pt(),weight,100,0,100);
    FillHist(p.prefix+region+p.hprefix+"dirap"+p.suffix+suf,dilepton.Rapidity(),weight,60,-3,3);
    if(p.bjets.size()){
      FillHist(p.prefix+region+p.hprefix+"b0pt"+p.suffix+suf,p.bjets.at(0).Pt(),weight,100,0,200);
      FillHist(p.prefix+region+p.hprefix+"b0eta"+p.suffix+suf,p.bjets.at(0).Eta(),weight,60,-3,3);
      FillHist(p.prefix+region+p.hprefix+"drbl"+p.suffix+suf,TMath::Min(p.bjets.at(0).DeltaR(*p.lepton0),p.bjets.at(0).DeltaR(*p.lepton1)),weight,60,-3,3);      
    }
    if(p.jets.size()){
      FillHist(p.prefix+region+p.hprefix+"j0pt"+p.suffix+suf,p.jets.at(0).Pt(),weight,100,0,200);
      FillHist(p.prefix+region+p.hprefix+"j0eta"+p.suffix+suf,p.jets.at(0).Eta(),weight,60,-3,3);
      FillHist(p.prefix+region+p.hprefix+"drjl"+p.suffix+suf,TMath::Min(p.jets.at(0).DeltaR(*p.lepton0),p.jets.at(0).DeltaR(*p.lepton1)),weight,60,-3,3);
    }      
    for(int i=0,n=p.leptons.size();i<n;i++){
      double eta=fabs(p.leptons.at(i)->Eta());
      double pt=p.leptons.at(i)->Pt();
      double riso=p.leptons.at(i)->RelIso();
      double energy=p.leptons.at(i)->E();
      FillHist(p.prefix+region+p.hprefix+Form("l%detapt",i)+p.suffix+suf,eta,pt,weight,netabin,etabins,nptbin,ptbins);
      FillHist(p.prefix+region+p.hprefix+Form("l%detae",i)+p.suffix+suf,eta,energy,weight,netabin,etabins,nptbin,ptbins);
      FillHist(p.prefix+region+p.hprefix+Form("l%driso",i)+p.suffix+suf,riso,weight,100,0,2);
      FillHist(p.prefix+region+p.hprefix+"lpt"+p.suffix+suf,pt,weight,nptbin,ptbins);
      FillHist(p.prefix+region+p.hprefix+Form("l%dpt",i)+p.suffix+suf,pt,weight,nptbin,ptbins);
      FillHist(p.prefix+region+p.hprefix+"leta"+p.suffix+suf,eta,weight,netabin,etabins);
      FillHist(p.prefix+region+p.hprefix+Form("l%deta",i)+p.suffix+suf,eta,weight,netabin,etabins);
    }
    if((p.channel=="EE"||p.channel=="MM")&&p.prefix.Contains("cpt/")){
      TString prefix=p.prefix+region;
      prefix.ReplaceAll("EE","ee");
      prefix.ReplaceAll("MM","mm");
      prefix.ReplaceAll("cpt/","");
      map<TString,int> map_fakesuf={{"",0},{"_faketf_up",1},{"_faketf_down",-1}};
      for(auto [fakesuf,sys]:map_fakesuf){
	double tf=GetFakeTF(p,"cpt "+region,sys);
	FillHist(prefix+"fakecpt_"+p.hprefix+"dimass"+p.suffix+suf+fakesuf,(*p.lepton0+*p.lepton1).M(),weight*tf,98,52,150);
	FillHist(prefix+"fakecpt_"+p.hprefix+"dimass_wide"+p.suffix+suf+fakesuf,(*p.lepton0+*p.lepton1).M(),weight*tf,nmbin,mbins);
	FillHist(prefix+"fakecpt_"+p.hprefix+"dipt"+p.suffix+suf+fakesuf,(*p.lepton0+*p.lepton1).Pt(),weight*tf,100,0,100);
	FillHist(prefix+"fakecpt_"+p.hprefix+"dirap"+p.suffix+suf+fakesuf,(*p.lepton0+*p.lepton1).Rapidity(),weight*tf,60,-3,3);
	FillHist(prefix+"fakecpt_"+p.hprefix+"lpt"+p.suffix+suf+fakesuf,p.lepton0->Pt(),weight*tf,nptbin,ptbins);
	FillHist(prefix+"fakecpt_"+p.hprefix+"lpt"+p.suffix+suf+fakesuf,p.lepton1->Pt(),weight*tf,nptbin,ptbins);
	FillHist(prefix+"fakecpt_"+p.hprefix+"l0pt"+p.suffix+suf+fakesuf,p.lepton0->Pt(),weight*tf,nptbin,ptbins);
	FillHist(prefix+"fakecpt_"+p.hprefix+"l1pt"+p.suffix+suf+fakesuf,p.lepton1->Pt(),weight*tf,nptbin,ptbins);
	FillHist(prefix+"fakecpt_"+p.hprefix+"leta"+p.suffix+suf+fakesuf,fabs(p.lepton0->Eta()),weight*tf,netabin,etabins);
	FillHist(prefix+"fakecpt_"+p.hprefix+"leta"+p.suffix+suf+fakesuf,fabs(p.lepton1->Eta()),weight*tf,netabin,etabins);
	FillHist(prefix+"fakecpt_"+p.hprefix+"l0eta"+p.suffix+suf+fakesuf,fabs(p.lepton0->Eta()),weight*tf,netabin,etabins);
	FillHist(prefix+"fakecpt_"+p.hprefix+"l1eta"+p.suffix+suf+fakesuf,fabs(p.lepton1->Eta()),weight*tf,netabin,etabins);
      }
    }else if((p.channel=="EE"||p.channel=="MM")&&p.prefix.Contains("mpt/")){
      TString prefix=p.prefix+region;
      prefix.ReplaceAll("EE","ee");
      prefix.ReplaceAll("MM","mm");
      prefix.ReplaceAll("mpt/","");
      map<TString,int> map_fakesuf={{"",0},{"_faketf_up",1},{"_faketf_down",-1}};
      for(auto [fakesuf,sys]:map_fakesuf){
	double tf=GetFakeTF(p,"mpt "+region,sys);
	FillHist(prefix+"fakempt_"+p.hprefix+"dimass"+p.suffix+suf+fakesuf,(*p.lepton0+*p.lepton1).M(),weight*tf,98,52,150);
	FillHist(prefix+"fakempt_"+p.hprefix+"dimass_wide"+p.suffix+suf+fakesuf,(*p.lepton0+*p.lepton1).M(),weight*tf,nmbin,mbins);
	FillHist(prefix+"fakempt_"+p.hprefix+"dipt"+p.suffix+suf+fakesuf,(*p.lepton0+*p.lepton1).Pt(),weight*tf,100,0,100);
	FillHist(prefix+"fakempt_"+p.hprefix+"dirap"+p.suffix+suf+fakesuf,(*p.lepton0+*p.lepton1).Rapidity(),weight*tf,60,-3,3);
	FillHist(prefix+"fakempt_"+p.hprefix+"lpt"+p.suffix+suf+fakesuf,p.lepton0->Pt(),weight*tf,nptbin,ptbins);
	FillHist(prefix+"fakempt_"+p.hprefix+"lpt"+p.suffix+suf+fakesuf,p.lepton1->Pt(),weight*tf,nptbin,ptbins);
	FillHist(prefix+"fakempt_"+p.hprefix+"l0pt"+p.suffix+suf+fakesuf,p.lepton0->Pt(),weight*tf,nptbin,ptbins);
	FillHist(prefix+"fakempt_"+p.hprefix+"l1pt"+p.suffix+suf+fakesuf,p.lepton1->Pt(),weight*tf,nptbin,ptbins);
	FillHist(prefix+"fakempt_"+p.hprefix+"leta"+p.suffix+suf+fakesuf,fabs(p.lepton0->Eta()),weight*tf,netabin,etabins);
	FillHist(prefix+"fakempt_"+p.hprefix+"leta"+p.suffix+suf+fakesuf,fabs(p.lepton1->Eta()),weight*tf,netabin,etabins);
	FillHist(prefix+"fakempt_"+p.hprefix+"l0eta"+p.suffix+suf+fakesuf,fabs(p.lepton0->Eta()),weight*tf,netabin,etabins);
	FillHist(prefix+"fakempt_"+p.hprefix+"l1eta"+p.suffix+suf+fakesuf,fabs(p.lepton1->Eta()),weight*tf,netabin,etabins);
      }
    }else if((p.channel=="EE"||p.channel=="MM")){
      TString prefix=p.prefix+region;
      prefix.ReplaceAll("EE","ee");
      prefix.ReplaceAll("MM","mm");
      map<TString,int> map_fakesuf={{"",0},{"_faketf_up",1},{"_faketf_down",-1}};
      for(auto [fakesuf,sys]:map_fakesuf){
	double tf=GetFakeTF(p,region,sys);
	FillHist(prefix+"fake_"+p.hprefix+"dimass"+p.suffix+suf+fakesuf,(*p.lepton0+*p.lepton1).M(),weight*tf,98,52,150);
	FillHist(prefix+"fake_"+p.hprefix+"dimass_wide"+p.suffix+suf+fakesuf,(*p.lepton0+*p.lepton1).M(),weight*tf,nmbin,mbins);
	FillHist(prefix+"fake_"+p.hprefix+"dipt"+p.suffix+suf+fakesuf,(*p.lepton0+*p.lepton1).Pt(),weight*tf,100,0,100);
	FillHist(prefix+"fake_"+p.hprefix+"dirap"+p.suffix+suf+fakesuf,(*p.lepton0+*p.lepton1).Rapidity(),weight*tf,60,-3,3);
	FillHist(prefix+"fake_"+p.hprefix+"lpt"+p.suffix+suf+fakesuf,p.lepton0->Pt(),weight*tf,nptbin,ptbins);
	FillHist(prefix+"fake_"+p.hprefix+"lpt"+p.suffix+suf+fakesuf,p.lepton1->Pt(),weight*tf,nptbin,ptbins);
	FillHist(prefix+"fake_"+p.hprefix+"l0pt"+p.suffix+suf+fakesuf,p.lepton0->Pt(),weight*tf,nptbin,ptbins);
	FillHist(prefix+"fake_"+p.hprefix+"l1pt"+p.suffix+suf+fakesuf,p.lepton1->Pt(),weight*tf,nptbin,ptbins);
	FillHist(prefix+"fake_"+p.hprefix+"leta"+p.suffix+suf+fakesuf,fabs(p.lepton0->Eta()),weight*tf,netabin,etabins);
	FillHist(prefix+"fake_"+p.hprefix+"leta"+p.suffix+suf+fakesuf,fabs(p.lepton1->Eta()),weight*tf,netabin,etabins);
	FillHist(prefix+"fake_"+p.hprefix+"l0eta"+p.suffix+suf+fakesuf,fabs(p.lepton0->Eta()),weight*tf,netabin,etabins);
	FillHist(prefix+"fake_"+p.hprefix+"l1eta"+p.suffix+suf+fakesuf,fabs(p.lepton1->Eta()),weight*tf,netabin,etabins);

	tf=GetFakeTF(p,"lj",sys);
	FillHist(prefix+"fakelj_"+p.hprefix+"dimass"+p.suffix+suf+fakesuf,(*p.lepton0+*p.lepton1).M(),weight*tf,98,52,150);
	FillHist(prefix+"fakelj_"+p.hprefix+"dimass_wide"+p.suffix+suf+fakesuf,(*p.lepton0+*p.lepton1).M(),weight*tf,nmbin,mbins);
	FillHist(prefix+"fakelj_"+p.hprefix+"dipt"+p.suffix+suf+fakesuf,(*p.lepton0+*p.lepton1).Pt(),weight*tf,100,0,100);
	FillHist(prefix+"fakelj_"+p.hprefix+"dirap"+p.suffix+suf+fakesuf,(*p.lepton0+*p.lepton1).Rapidity(),weight*tf,60,-3,3);
	FillHist(prefix+"fakelj_"+p.hprefix+"lpt"+p.suffix+suf+fakesuf,p.lepton0->Pt(),weight*tf,nptbin,ptbins);
	FillHist(prefix+"fakelj_"+p.hprefix+"lpt"+p.suffix+suf+fakesuf,p.lepton1->Pt(),weight*tf,nptbin,ptbins);
	FillHist(prefix+"fakelj_"+p.hprefix+"l0pt"+p.suffix+suf+fakesuf,p.lepton0->Pt(),weight*tf,nptbin,ptbins);
	FillHist(prefix+"fakelj_"+p.hprefix+"l1pt"+p.suffix+suf+fakesuf,p.lepton1->Pt(),weight*tf,nptbin,ptbins);
	FillHist(prefix+"fakelj_"+p.hprefix+"leta"+p.suffix+suf+fakesuf,fabs(p.lepton0->Eta()),weight*tf,netabin,etabins);
	FillHist(prefix+"fakelj_"+p.hprefix+"leta"+p.suffix+suf+fakesuf,fabs(p.lepton1->Eta()),weight*tf,netabin,etabins);
	FillHist(prefix+"fakelj_"+p.hprefix+"l0eta"+p.suffix+suf+fakesuf,fabs(p.lepton0->Eta()),weight*tf,netabin,etabins);
	FillHist(prefix+"fakelj_"+p.hprefix+"l1eta"+p.suffix+suf+fakesuf,fabs(p.lepton1->Eta()),weight*tf,netabin,etabins);
      }
    }
  }
}
vector<Muon> FakeAnalyzer::ToConePt(vector<Muon> muons){
  for(auto& muon:muons){
    muon*=1+muon.RelIso();
  }
  std::sort(muons.begin(),muons.end(),PtComparing);
  return muons;
}
vector<Electron> FakeAnalyzer::ToConePt(vector<Electron> electrons){
  for(auto& electron:electrons){
    electron*=1+electron.RelIso();
  }
  std::sort(electrons.begin(),electrons.end(),PtComparing);
  return electrons;
}
vector<Muon> FakeAnalyzer::ToModifiedPt(vector<Muon> muons){
  for(auto& muon:muons){
    muon*=1+muon.RelIso()*TMath::Max(0.,TMath::Min(1.,(muon.Pt()-30)/30));
  }
  std::sort(muons.begin(),muons.end(),PtComparing);
  return muons;
}
vector<Electron> FakeAnalyzer::ToModifiedPt(vector<Electron> electrons){
  for(auto& electron:electrons){
    electron*=1+electron.RelIso()*TMath::Max(0.,TMath::Min(1.,(electron.Pt()-30)/30));
  }
  std::sort(electrons.begin(),electrons.end(),PtComparing);
  return electrons;
}
SMPAnalyzerCore::Parameter FakeAnalyzer::MakeParameter(TString channel,TString option){
  Parameter p=SMPAnalyzerCore::MakeParameter(channel,option);
  if(IsDATA&&p.channel.Contains("j")){
    double maxlumi=0;
    for(TString trigger:p.triggers){
      if(_event.PassTrigger(trigger)){
	lumi=_event.GetTriggerLumi(trigger);
	if(lumi>maxlumi) maxlumi=lumi;
      }
    }
    double fulllumi=_event.GetTriggerLumi("Full");
    p.w.lumiweight*=fulllumi/maxlumi;
  }
  return p;
}
