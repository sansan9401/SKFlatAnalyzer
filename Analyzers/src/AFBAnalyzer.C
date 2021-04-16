#include "AFBAnalyzer.h"

void AFBAnalyzer::initializeAnalyzer(){
  SMPAnalyzerCore::initializeAnalyzer(); //setup zpt roc z0 
  //SetupToy(100);
  SetupCosThetaWeight();
  IsNominalRun=!HasFlag("SYS")&&!HasFlag("PDFSYS");
  
  if(HasFlag("bjet")||HasFlag("nobjet")){
    vector<JetTagging::Parameters> jtps={JetTagging::Parameters(JetTagging::DeepCSV,JetTagging::Medium,JetTagging::mujets,JetTagging::mujets),
                                         JetTagging::Parameters(JetTagging::DeepJet,JetTagging::Medium,JetTagging::mujets,JetTagging::mujets)};
    mcCorr->SetJetTaggingParameters(jtps);
  }
  if(fChain->GetListOfFiles()->GetEntries()){
    TString filename=fChain->GetListOfFiles()->At(0)->GetTitle();
    if(filename.Contains("SkimTree_")) IsSkimmed=true;
    else IsSkimmed=false;
  }else{
    cout<<"[AFBAnalyzer::initializeAnalyzer] no input file"<<endl;
    exit(EXIT_FAILURE);
  }
  if(!IsSkimmed&&!HasFlag("ALL")){
    fChain->SetBranchStatus("pfMET_*",false);
    fChain->SetBranchStatus("jet_*",false);
    fChain->SetBranchStatus("fatjet_*",false);
    fChain->SetBranchStatus("electron_*",false);
    fChain->SetBranchStatus("muon_*",false);
    fChain->SetBranchStatus("photon_*",false);
  }
}
void AFBAnalyzer::executeEvent(){
  //GetToyWeight();
  costhetaweight=1.;
  costhetaweight_up=1.;
  costhetaweight_down=1.;
  if(IsDYSample){
    //////////////////////// Check LHE /////////////////////////
    if(abs(lhe_l0.ID())!=15){
      Parameter p;
      double letacut=2.4;
      if(abs(lhe_l0.ID())==11){
	p=MakeParameter("ee");
	p.c.lepton0pt=25;
	p.c.lepton1pt=15;
      }else if(abs(lhe_l0.ID())==13){
	p=MakeParameter("mm");
	p.c.lepton0pt=20;
	p.c.lepton1pt=10;
      }else{
	cout<<"[AFBAnalyzer::executeEvent()] something is wrong l0.ID="<<abs(lhe_l0.ID())<<endl;
	vector<LHE> lhes=GetLHEs();
	for(auto& lhe:lhes) lhe.Print();
	exit(EXIT_FAILURE);
      }
      
      //////////////////////// GEN /////////////////////////
      TLorentzVector gen_Z=gen_l0+gen_l1;
      double gen_Zmass=gen_Z.M();
      double gen_Zrap=gen_Z.Rapidity();
      double gen_Zpt=gen_Z.Pt();
      double gen_cost_correct=-999;
      if(gen_p0.PID()==21){
	if(gen_p1.PID()==21) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,0);
	else if(gen_p1.PID()>0) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,-1);
	else if(gen_p1.PID()<0) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,1);
      }else if(gen_p0.PID()>0){
	if(gen_p1.PID()==21) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,1);
	else if(gen_p1.PID()>0) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,0);
	else if(gen_p1.PID()<0) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,1);
      }else if(gen_p0.PID()<0){
	if(gen_p1.PID()==21) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,-1);
	else if(gen_p1.PID()>0) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,-1);
	else if(gen_p1.PID()<0) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,0);
      }
      if(gen_cost_correct==-999){
	cout<<"wrong pid for parton"<<endl;
	exit(EXIT_FAILURE);
      }
      costhetaweight=GetCosThetaWeight(gen_Zmass,gen_Zpt,gen_cost_correct,"_pdg");
      costhetaweight_up=GetCosThetaWeight(gen_Zmass,gen_Zpt,gen_cost_correct,"_up");
      costhetaweight_down=GetCosThetaWeight(gen_Zmass,gen_Zpt,gen_cost_correct,"_down");
      
      map<TString,double> map_weight;
      map_weight[""]=p.w.lumiweight*p.w.zptweight*costhetaweight;
      map_weight["_noweight"]=p.w.lumiweight;
      map_weight["_nozptweight"]=p.w.lumiweight*costhetaweight;
      map_weight["_nocosthetaweight"]=p.w.lumiweight*p.w.zptweight;

      //////////////// Fill LHE,Gen hists //////////////////////
      if(IsNominalRun&&!IsSkimmed&&!HasFlag("ALL")){
	FillHistsAFB(p.prefix,"lhe_","",(Particle*)&lhe_l0,(Particle*)&lhe_l1,map_weight);
	FillHistsAFB(p.prefix,"gen_","",(Particle*)&gen_l0,(Particle*)&gen_l1,map_weight);
	FillHistsAFB(p.prefix,"gen_","_dressed",(Particle*)&gen_l0_dressed,(Particle*)&gen_l1_dressed,map_weight);
	FillHistsAFB(p.prefix,"gen_","_bare",(Particle*)&gen_l0_bare,(Particle*)&gen_l1_bare,map_weight);
	if(gen_l0.Pt()>p.c.lepton0pt&&gen_l1.Pt()>p.c.lepton1pt&&fabs(gen_l0.Eta())<letacut&&fabs(gen_l1.Eta())<letacut){
	  FillHistsAFB(p.prefix,"genfid_","",(Particle*)&gen_l0,(Particle*)&gen_l1,map_weight);
	}
	if(gen_l0_dressed.Pt()>p.c.lepton0pt&&gen_l1_dressed.Pt()>p.c.lepton1pt&&fabs(gen_l0_dressed.Eta())<letacut&&fabs(gen_l1_dressed.Eta())<letacut){
	  FillHistsAFB(p.prefix,"genfid_","_dressed",(Particle*)&gen_l0_dressed,(Particle*)&gen_l1_dressed,map_weight);
	}
	if(gen_l0_bare.Pt()>p.c.lepton0pt&&gen_l1_bare.Pt()>p.c.lepton1pt&&fabs(gen_l0_bare.Eta())<letacut&&fabs(gen_l1_bare.Eta())<letacut){
	  FillHistsAFB(p.prefix,"genfid_","_bare",(Particle*)&gen_l0_bare,(Particle*)&gen_l1_bare,map_weight);
	}
	FillHist(p.prefix+"gen_costhetaCS_correct",gen_Zmass,gen_Zrap,gen_Zpt,gen_cost_correct,map_weight,afb_mbinnum,(double*)afb_mbin,afb_ybinnum,(double*)afb_ybin,afb_ptbinnum,(double*)afb_ptbin,20,-1,1);
      }
    }
  }

  if(!IsSkimmed&&!HasFlag("ALL")) return;

  ///////////////// RECO level /////////////////////
  if(!IsDATA||DataStream.Contains("SingleMuon")){
    executeEventWithParameter(MakeParameter("me"));
    executeEventWithParameter(MakeParameter("mM"));
    executeEventWithParameter(MakeParameter("mu"));
  }
  if(!IsDATA||DataStream.Contains("DoubleMuon")){
    executeEventWithParameter(MakeParameter("mm"));
    executeEventWithParameter(MakeParameter("MM"));
  }
  if(!IsDATA||DataStream.Contains("SingleElectron")||DataStream.Contains("EGamma")){
    executeEventWithParameter(MakeParameter("em"));
    executeEventWithParameter(MakeParameter("eE"));
    executeEventWithParameter(MakeParameter("el"));
  }
  if(!IsDATA||DataStream.Contains("DoubleEG")||DataStream.Contains("EGamma")){
    executeEventWithParameter(MakeParameter("ee"));
    executeEventWithParameter(MakeParameter("EE"));
  }
}
SMPAnalyzerCore::Parameter AFBAnalyzer::MakeParameter(TString key){
  Parameter p=SMPAnalyzerCore::MakeParameter(key);

  if(IsNominalRun) p.weightbit|=NominalWeight;
  if(HasFlag("SYS")&&!IsDATA) p.weightbit|=SystematicWeight;
  if(HasFlag("PDFSYS")&&!IsDATA) p.weightbit|=PDFWeight;

  if(HasFlag("bjet")) p.prefix+="bjet/";
  else if(HasFlag("nobjet")) p.prefix+="nobjet/";
  if(HasFlag("highmet")) p.prefix+="highmet/";

  return p;
}
bool AFBAnalyzer::PassSelection(Parameter& p){
  if(p.prefix.Contains("highmet")){
    if(pfMET_Type1_pt<60) return false;
    if(IsNominalRun) FillCutflow(p.prefix+p.hprefix+"cutflow","METCut",p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight);
  }

  int n_bjet=0;
  if(p.prefix.Contains("bjet")||p.prefix.Contains("nobjet")){
    std::vector<Jet> jets=GetJets("tightLepVeto",30,2.7);
    std::sort(jets.begin(),jets.end(),PtComparing);

    JetTagging::Parameters jtp = JetTagging::Parameters(JetTagging::DeepCSV,JetTagging::Medium,JetTagging::mujets,JetTagging::mujets);
    for(const auto& jet:jets)
      if(mcCorr->IsBTagged_2a(jtp,jet))
	n_bjet++;
    
    if(p.prefix.Contains("bjet")&&!n_bjet) return false;
    if(p.prefix.Contains("nobjet")&&n_bjet) return false;
    if(IsNominalRun) FillCutflow(p.prefix+p.hprefix+"cutflow","BJetCut",p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight);
  }

  if(!SMPAnalyzerCore::PassSelection(p)) return false;  
  return true;
}

void AFBAnalyzer::FillHists(Parameter& p){
  ///////////////////////map_weight//////////////////
  map<TString,double> map_weight;
  if(p.weightbit&NominalWeight){
    map_weight[""]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF*p.w.CFSF;
    map_weight["_noCFSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF;
    map_weight["_CFSF_up"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF*p.w.CFSF_up;
    map_weight["_CFSF_down"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF*p.w.CFSF_down;
  }
  if(p.weightbit&SystematicWeight){
    map_weight["_noPUweight"]=p.w.lumiweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF;
    map_weight["_PUweight_up"]=p.w.lumiweight*p.w.PUweight_up*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF;
    map_weight["_PUweight_down"]=p.w.lumiweight*p.w.PUweight_down*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF;
    
    map_weight["_noprefireweight"]=p.w.lumiweight*p.w.PUweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF;
    map_weight["_prefireweight_up"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight_up*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF;
    map_weight["_prefireweight_down"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight_down*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF;
    
    map_weight["_nozptweight"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF;
    
    map_weight["_noz0weight"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF;

    map_weight["_nocosthetaweight"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF;
    map_weight["_costhetaweight_up"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight_up*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF;
    map_weight["_costhetaweight_down"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight_down*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF;
    
    map_weight["_noefficiencySF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight;
    
    map_weight["_noRECOSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.IDSF*p.w.ISOSF*p.w.triggerSF;
    map_weight["_RECOSF_up"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF_up*p.w.IDSF*p.w.ISOSF*p.w.triggerSF;
    map_weight["_RECOSF_down"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF_down*p.w.IDSF*p.w.ISOSF*p.w.triggerSF;
    
    map_weight["_noIDSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.ISOSF*p.w.triggerSF;
    map_weight["_IDSF_up"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF_up*p.w.ISOSF*p.w.triggerSF;
    map_weight["_IDSF_down"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF_down*p.w.ISOSF*p.w.triggerSF;
    
    map_weight["_noISOSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.triggerSF;
    map_weight["_ISOSF_up"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF_up*p.w.triggerSF;
    map_weight["_ISOSF_down"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF_down*p.w.triggerSF;
    
    map_weight["_notriggerSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF;
    map_weight["_triggerSF_up"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF_up;
    map_weight["_triggerSF_down"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF_down;
    
  }
  if(p.weightbit&PDFWeight){
    for(unsigned int i=0;i<PDFWeights_Scale->size();i++){
      map_weight[Form("_scalevariation%d",i)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF*PDFWeights_Scale->at(i);
    }
    for(unsigned int i=0;i<PDFWeights_Error->size();i++){
      map_weight[Form("_pdf%d",i)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF*PDFWeights_Error->at(i);
    }
    if(PDFWeights_AlphaS->size()==2){
      map_weight["_alphaS_up"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF*PDFWeights_AlphaS->at(0);
      map_weight["_alphaS_down"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF*PDFWeights_AlphaS->at(1);
    }
  }
  
  ///////////////////////fill hists///////////////////////
  if(HasFlag("TOY")) FillHistsToy(p.prefix,p.hprefix,p.suffix,(Particle*)p.lepton0,(Particle*)p.lepton1,map_weight);
  else{
    FillHistsAFB(p.prefix,p.hprefix,p.suffix,(Particle*)p.lepton0,(Particle*)p.lepton1,map_weight);
    if(IsDYSample&&p.hprefix==""&&IsNominalRun){
      vector<Gen> gens=GetGens();
      Gen truth_l0=GetGenMatchedLepton(*p.lepton0,gens);
      Gen truth_l1=GetGenMatchedLepton(*p.lepton1,gens);
      if(!truth_l0.IsEmpty()&&!truth_l1.IsEmpty()) 
	FillHistsAFB(p.prefix,"truth_",p.suffix,(Particle*)&truth_l0,(Particle*)&truth_l1,map_weight);
      //else cout<<"no matching"<<endl;
    }
  }
}

AFBAnalyzer::AFBAnalyzer(){}
AFBAnalyzer::~AFBAnalyzer(){
  //DeleteToy();
  DeleteCosThetaWeight();
}
double AFBAnalyzer::GetCosThetaCS(const Particle *p0,const Particle *p1,int direction){
  const TLorentzVector *l0,*l1;
  if(p0->Charge()<0&&p1->Charge()>0){
    l0=p0;
    l1=p1;
  }else if(p0->Charge()>0&&p1->Charge()<0){
    l0=p1;
    l1=p0;
  }else if(strcmp(p0->ClassName(),"LHE")==0){ 
    if(((LHE*)p0)->ID()>0&&((LHE*)p1)->ID()<0){
      l0=p0;
      l1=p1;
    }else if(((LHE*)p0)->ID()<0&&((LHE*)p1)->ID()>0){
      l0=p1;
      l1=p0;
    }else{
      if(gRandom->Rndm()<0.5){
	l0=p0;
	l1=p1;
      }else{
	l0=p1;
	l1=p0;
      }      
    } 
  }else{
    if(gRandom->Rndm()<0.5){
      l0=p0;
      l1=p1;
    }else{
      l0=p1;
      l1=p0;
    }      
  }

  TLorentzVector dilepton=*l0+*l1;
  double l0pp=(l0->E()+l0->Pz())/sqrt(2);
  double l0pm=(l0->E()-l0->Pz())/sqrt(2);
  double l1pp=(l1->E()+l1->Pz())/sqrt(2);
  double l1pm=(l1->E()-l1->Pz())/sqrt(2);
  double dimass=dilepton.M();
  double dipt=dilepton.Pt();
  if(direction==0) direction=dilepton.Pz()>0?1:-1;
  return direction*2*(l0pp*l1pm-l0pm*l1pp)/sqrt(dimass*dimass*(dimass*dimass+dipt*dipt));
} 
double AFBAnalyzer::GetCosTheta(const vector<Lepton*>& leps,const vector<Jet>& jets,TString option,double fcut){
  TLorentzVector *l0,*l1;
  if(leps.at(0)->Charge()<0&&leps.at(1)->Charge()>0){
    l0=leps.at(0);
    l1=leps.at(1);
  }else if(leps.at(0)->Charge()>0&&leps.at(1)->Charge()<0){
    l0=leps.at(1);
    l1=leps.at(0);
  }else{
    if(gRandom->Rndm()<0.5){
      l0=leps.at(0);
      l1=leps.at(1);
    }else{
      l0=leps.at(1);
      l1=leps.at(0);
    }      
  }
  TLorentzVector dilepton=*l0+*l1;
  TLorentzVector j0;
  if(jets.size()>0){
    if((jets[0]+dilepton).Pt()/dilepton.Pt()<fcut)
      j0=jets[0];
  }
  TLorentzVector system=dilepton+j0;
  double direction=dilepton.Pz()/fabs(dilepton.Pz());
  double systemmass=system.M();
  //double systempt=system.Pt();
  double systempz=system.Pz();
  TLorentzVector p0(0,0,(sqrt(systemmass*systemmass+systempz*systempz)+systempz)/2.,(sqrt(systemmass*systemmass+systempz*systempz)+systempz)/2.);
  TLorentzVector p1(0,0,-(sqrt(systemmass*systemmass+systempz*systempz)-systempz)/2.,(sqrt(systemmass*systemmass+systempz*systempz)-systempz)/2.);
  TLorentzVector particle=*l0;
  TVector3 b1=system.BoostVector();
  dilepton.Boost(-b1);particle.Boost(-b1);
  p0.Boost(-b1);p1.Boost(-b1);j0.Boost(-b1);
  if(j0.E()>0){
    if(p0.Angle(j0.Vect())<p1.Angle(j0.Vect())){
      if(option.Contains("C")) p0-=j0;
      if(option.Contains("B")) direction=-1;
    }else{
      if(option.Contains("C")) p1-=j0;
      if(option.Contains("B")) direction=1;
    }
  }
  TVector3 b2=dilepton.BoostVector();
  p0.Boost(-b2);p1.Boost(-b2);particle.Boost(-b2);
  return direction*cos(particle.Angle(p0.Vect().Unit()-p1.Vect().Unit()));
}
void AFBAnalyzer::FillHistsToy(TString pre,TString hpre,TString suf,Particle* l0,Particle* l1,map<TString,double> map_weight){
  int n_toy=toy_random.size();
  for(int i=0;i<n_toy;i++) FillHistsAFB(pre,hpre,suf+Form("_toy%d",i),l0,l1,Multiply(map_weight,toy_weight[i]));
}
void AFBAnalyzer::FillHistsAFB(TString pre,TString hpre,TString suf,Particle* l0,Particle* l1,map<TString,double> map_weight){
  TLorentzVector dilepton=(*l0)+(*l1);
  double dimass=dilepton.M();
  double dirap=dilepton.Rapidity();
  double dipt=dilepton.Pt();

  double cost=GetCosThetaCS(l0,l1);
  double h=0.5*pow(dipt/dimass,2)/(1+pow(dipt/dimass,2))*(1-3*cost*cost);
  double den_weight=0.5*fabs(cost)/pow(1+cost*cost+h,2);
  double num_weight=0.5*cost*cost/pow(1+cost*cost+h,3);
  FillHist(pre+hpre+"costhetaCS"+suf,dimass,dirap,dipt,cost,map_weight,afb_mbinnum,(double*)afb_mbin,afb_ybinnum,(double*)afb_ybin,afb_ptbinnum,(double*)afb_ptbin,20,-1,1);
  //FillHist(pre+hpre+"costhetaCS_den"+suf,dimass,dirap,dipt,cost,Multiply(map_weight,den_weight),afb_mbinnum,(double*)afb_mbin,afb_ybinnum,(double*)afb_ybin,afb_ptbinnum,(double*)afb_ptbin,20,-1,1);
  //FillHist(pre+hpre+"costhetaCS_num"+suf,dimass,dirap,dipt,cost,Multiply(map_weight,num_weight),afb_mbinnum,(double*)afb_mbin,afb_ybinnum,(double*)afb_ybin,afb_ptbinnum,(double*)afb_ptbin,20,-1,1);

  if(!HasFlag("PDFSYS")){
    FillHist(pre+hpre+"l0pt"+suf,dimass,dirap,dipt,l0->Pt(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,lptbinnum,(double*)lptbin);
    FillHist(pre+hpre+"l1pt"+suf,dimass,dirap,dipt,l1->Pt(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,lptbinnum,(double*)lptbin);
    FillHist(pre+hpre+"lpt"+suf,dimass,dirap,dipt,l0->Pt(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,lptbinnum,(double*)lptbin);
    FillHist(pre+hpre+"lpt"+suf,dimass,dirap,dipt,l1->Pt(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,lptbinnum,(double*)lptbin);
    
    FillHist(pre+hpre+"l0eta"+suf,dimass,dirap,dipt,l0->Eta(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,60,-3,3);
    FillHist(pre+hpre+"l1eta"+suf,dimass,dirap,dipt,l1->Eta(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,60,-3,3);
    FillHist(pre+hpre+"leta"+suf,dimass,dirap,dipt,l0->Eta(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,60,-3,3);
    FillHist(pre+hpre+"leta"+suf,dimass,dirap,dipt,l1->Eta(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,60,-3,3);

    if(!hpre.Contains("gen")&&!hpre.Contains("lhe")&&!hpre.Contains("truth")){
      FillHist(pre+hpre+"z0"+suf,dimass,dirap,dipt,vertex_Z,map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,120,-15,15);
      FillHist(pre+hpre+"met"+suf,dimass,dirap,dipt,pfMET_Type1_pt,map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,100,0,200);  
      FillHist(pre+hpre+"nPV"+suf,dimass,dirap,dipt,nPV,map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,60,0,60);  
      FillHist(pre+hpre+"nPileUp"+suf,dimass,dirap,dipt,nPileUp,map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,60,0,60);  
    }
  }
}
void AFBAnalyzer::FillHardHists(TString pre,TString suf,const Gen& genparton0,const Gen& genparton1,const Gen& genhardl0,const Gen& genhardl1,const Gen& genhardj0,double w){
  Gen genhardl=genhardl0.PID()>0?genhardl0:genhardl1;
  TLorentzVector genZ=genhardl0+genhardl1;
  TLorentzVector genpp=genparton0+genparton1;

  FillHist(pre+"cos_l_p0"+suf,cos(genhardl.Angle(genparton0.Vect())),w,100,-1,1);
  FillHist(pre+"cos_l_p0_asym"+suf,cos(genhardl.Angle(genparton0.Vect())),w/2,100,-1,1);
  FillHist(pre+"cos_l_p0_asym"+suf,-1.*cos(genhardl.Angle(genparton0.Vect())),-w/2,100,-1,1);
  if(cos(genhardl.Angle(genparton0.Vect()))>0) FillHist(pre+"cos_l_p0_forward"+suf,genZ.M(),w,fine_mbinnum,(double*)fine_mbin);
  else FillHist(pre+"cos_l_p0_backward"+suf,genZ.M(),w,fine_mbinnum,(double*)fine_mbin);

  FillHist(pre+"cos_l_p1"+suf,cos(genhardl.Angle(genparton1.Vect())),w,100,-1,1);
  FillHist(pre+"cos_l_p1_asym"+suf,cos(genhardl.Angle(genparton1.Vect())),w/2,100,-1,1);
  FillHist(pre+"cos_l_p1_asym"+suf,-1.*cos(genhardl.Angle(genparton1.Vect())),-w/2,100,-1,1);
  if(cos(genhardl.Angle(genparton1.Vect()))>0) FillHist(pre+"cos_l_p1_forward"+suf,genZ.M(),w,fine_mbinnum,(double*)fine_mbin);
  else FillHist(pre+"cos_l_p1_backward"+suf,genZ.M(),w,fine_mbinnum,(double*)fine_mbin);
  
  FillHist(pre+"cos_Z_p0"+suf,cos(genZ.Angle(genparton0.Vect())),w,100,-1,1);
  FillHist(pre+"cos_Z_p0_asym"+suf,cos(genZ.Angle(genparton0.Vect())),w/2,100,-1,1);
  FillHist(pre+"cos_Z_p0_asym"+suf,-1.*cos(genZ.Angle(genparton0.Vect())),-w/2,100,-1,1);
  
  FillHist(pre+"cos_Z_p1"+suf,cos(genZ.Angle(genparton1.Vect())),w,100,-1,1);
  FillHist(pre+"cos_Z_p1_asym"+suf,cos(genZ.Angle(genparton1.Vect())),w/2,100,-1,1);
  FillHist(pre+"cos_Z_p1_asym"+suf,-1.*cos(genZ.Angle(genparton1.Vect())),-w/2,100,-1,1);
  
  
  FillHist(pre+"Zrap"+suf,genZ.Rapidity(),w,100,-5,5);
  FillHist(pre+"Zrap_asym"+suf,genZ.Rapidity(),w/2,100,-5,5);
  FillHist(pre+"Zrap_asym"+suf,-1.*genZ.Rapidity(),-w/2,100,-5,5);
  
  FillHist(pre+"pprap"+suf,genpp.Rapidity(),w,100,-5,5);
  FillHist(pre+"pprap_asym"+suf,genpp.Rapidity(),w/2,100,-5,5);
  FillHist(pre+"pprap_asym"+suf,-1.*genpp.Rapidity(),-w/2,100,-5,5);
  
  if(!genhardj0.IsEmpty()){
    FillHist(pre+"cos_l_j0"+suf,cos(genhardl.Angle(genhardj0.Vect())),w,100,-1,1);
    FillHist(pre+"cos_l_j0_asym"+suf,cos(genhardl.Angle(genhardj0.Vect())),w/2,100,-1,1);
    FillHist(pre+"cos_l_j0_asym"+suf,-1.*cos(genhardl.Angle(genhardj0.Vect())),-w/2,100,-1,1);
    if(cos(genhardl.Angle(genhardj0.Vect()))>0) FillHist(pre+"cos_l_j0_forward"+suf,genZ.M(),w,fine_mbinnum,(double*)fine_mbin);
    else FillHist(pre+"cos_l_j0_backward"+suf,genZ.M(),w,fine_mbinnum,(double*)fine_mbin);

    FillHist(pre+"cos_Z_j0"+suf,cos(genZ.Angle(genhardj0.Vect())),w,100,-1,1);
    FillHist(pre+"cos_Z_j0_asym"+suf,cos(genZ.Angle(genhardj0.Vect())),w/2,100,-1,1);
    FillHist(pre+"cos_Z_j0_asym"+suf,-1.*cos(genZ.Angle(genhardj0.Vect())),-w/2,100,-1,1);
    
    FillHist(pre+"cos_j0_p0"+suf,cos(genhardj0.Angle(genparton0.Vect())),w,100,-1,1);
    FillHist(pre+"cos_j0_p0_asym"+suf,cos(genhardj0.Angle(genparton0.Vect())),w/2,100,-1,1);
    FillHist(pre+"cos_j0_p0_asym"+suf,-1.*cos(genhardj0.Angle(genparton0.Vect())),-w/2,100,-1,1);
    
    FillHist(pre+"cos_j0_p1"+suf,cos(genhardj0.Angle(genparton1.Vect())),w,100,-1,1);
    FillHist(pre+"cos_j0_p1_asym"+suf,cos(genhardj0.Angle(genparton1.Vect())),w/2,100,-1,1);
    FillHist(pre+"cos_j0_p1_asym"+suf,-1.*cos(genhardj0.Angle(genparton1.Vect())),-w/2,100,-1,1);
    
    FillHist(pre+"jeta"+suf,genhardj0.Eta(),w,100,-5,5);
    FillHist(pre+"jeta_asym"+suf,genhardj0.Eta(),w/2,100,-5,5);
    FillHist(pre+"jeta_asym"+suf,-1.*genhardj0.Eta(),-w/2,100,-5,5);
  }
}
void AFBAnalyzer::SetupToy(int n_toy){
  DeleteToy();
  for(int i=0;i<n_toy;i++){
    toy_random.push_back(new TRandom3);
    toy_weight.push_back(-9999.);
  } 
}
void AFBAnalyzer::DeleteToy(){
  int n_toy=toy_random.size();
  for(int i=0;i<n_toy;i++) delete toy_random.at(i);
  toy_random.clear();
  toy_weight.clear();
}
void AFBAnalyzer::GetToyWeight(){
  int n_toy=toy_random.size();
  if(fChain->GetTree()->GetReadEntry()==0){
    TString filename=fChain->GetFile()->GetName();
    for(int i=0;i<n_toy;i++) toy_random[i]->SetSeed((filename+Form("%d",i)).MD5().Hash());
  }
  for(int i=0;i<n_toy;i++)
    toy_weight[i]=toy_random[i]->PoissonD(1.);
}

void AFBAnalyzer::FillHistToy(TString histname, double value, double weight, int n_bin, double x_min, double x_max){
  int n_toy=toy_random.size();
  for(int i=0;i<n_toy;i++) FillHist(histname+Form("_toy%d",i),value,weight*toy_weight[i],n_bin,x_min,x_max);
}
void AFBAnalyzer::FillHistToy(TString histname, double value, map<TString,double> weights, int n_bin, double x_min, double x_max){
  for(const auto& [suffix,weight]:weights) FillHistToy(histname+suffix,value,weight,n_bin,x_min,x_max);
}
void AFBAnalyzer::FillHistToy(TString histname, double value, double weight, int n_bin, double *xbins){ 
  int n_toy=toy_random.size();
  for(int i=0;i<n_toy;i++) FillHist(histname+Form("_toy%d",i),value,weight*toy_weight[i],n_bin,xbins);
}
void AFBAnalyzer::FillHistToy(TString histname, double value, map<TString,double> weights, int n_bin, double *xbins){ 
  for(const auto& [suffix,weight]:weights) FillHistToy(histname+suffix,value,weight,n_bin,xbins);
}
void AFBAnalyzer::FillHistToy(TString histname, double value_x, double value_y, double value_z, double weight, int n_binx, double *xbins, int n_biny, double *ybins, int n_binz, double *zbins){
  int n_toy=toy_random.size();
  for(int i=0;i<n_toy;i++) FillHist(histname+Form("_toy%d",i),value_x,value_y,value_z,weight,n_binx,xbins,n_biny,ybins,n_binz,zbins);
}
void AFBAnalyzer::FillHistToy(TString histname, double value_x, double value_y, double value_z, map<TString,double> weights, int n_binx, double *xbins, int n_biny, double *ybins, int n_binz, double *zbins){
  for(const auto& [suffix,weight]:weights) FillHistToy(histname+suffix,value_x,value_y,value_z,weight,n_binx,xbins,n_biny,ybins,n_binz,zbins);
}
void AFBAnalyzer::SetupCosThetaWeight(){
  cout<<"[AFBAnalyzer::SetupCosThetaWeight] Setup"<<endl;
  TString datapath=getenv("DATA_DIR");
  ifstream file_check(datapath+"/"+GetEra()+"/SMP/CosThetaWeight.root");
  bool isexist=file_check.is_open();
  file_check.close();
  if(!isexist){
    cout<<"[AFBAnalyzer::SetupCosThetaWeight] no CosThetaWeight.root"<<endl;
    return;
  }
  TFile fcost(datapath+"/"+GetEra()+"/SMP/CosThetaWeight.root");
  for(const auto&& key:*(fcost.GetListOfKeys())){
    TObject* obj=((TKey*)key)->ReadObj();
    if(!obj->InheritsFrom("TH3D")) continue;
    TH3D* hist=(TH3D*)obj;
    cout<<"[AFBAnalyzer::SetupCosThetaWeight] get "<<hist->GetName()<<endl;
    map_hist_cost[hist->GetName()]=hist;
    hist->SetDirectory(0);
  }
}
void AFBAnalyzer::DeleteCosThetaWeight(){
  for(auto& iter:map_hist_cost)
    if(iter.second) delete iter.second;
}
double AFBAnalyzer::GetCosThetaWeight(double mass,double pt,double cost,TString suffix){
  double val=1.;
  if(!IsDYSample) return val;
  TString MCName=MCSample;
  if(MCName.Contains(TRegexp("^DY[0-9]Jets$"))) MCName="DYJets";
  if(MCName.Contains(TRegexp("^DYJets_Pt-[0-9]*To[0-9Inf]*$"))) MCName="DYJets";
  if(MCName.Contains(TRegexp("^DYJets_M-[0-9]*to[0-9Inf]*$"))) MCName="DYJets";
  TString hname=MCName+suffix;
  auto it=map_hist_cost.find(hname);
  if(it!=map_hist_cost.end())
    val*=GetBinContentUser(it->second,mass,pt,cost,0);
  if(val==0) val=1.;
  return val;
}
