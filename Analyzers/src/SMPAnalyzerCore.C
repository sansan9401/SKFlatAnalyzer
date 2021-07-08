#include "SMPAnalyzerCore.h"

SMPAnalyzerCore::SMPAnalyzerCore(){}
SMPAnalyzerCore::~SMPAnalyzerCore(){
  for(auto& iter:map_hist_zpt){
    if(iter.second) delete iter.second;
  }
  if(roc) delete roc;
  if(rocele) delete rocele;
  if(hz0_data) delete hz0_data;
  if(hz0_mc) delete hz0_mc;
  for(std::map< TString, TH4D* >::iterator mapit = maphist_TH4D.begin(); mapit!=maphist_TH4D.end(); mapit++){
    delete mapit->second;
  }
  maphist_TH4D.clear();
  DeleteEfficiency();
  DeleteCFRate();
}

void SMPAnalyzerCore::initializeAnalyzer(){
  if(MaxEvent>0) reductionweight=1.*fChain->GetEntries()/MaxEvent;
  else reductionweight=1.;
  SetupEfficiency();
  SetupZptWeight();
  SetupRoccoR();
  SetupZ0Weight();
  SetupCFRate();
  IsDYSample=false;
  if(MCSample.Contains("DYJets")||MCSample.Contains("ZToEE")||MCSample.Contains("ZToMuMu")||MCSample.Contains(TRegexp("DY[0-9]Jets"))) IsDYSample=true;
}
void SMPAnalyzerCore::beginEvent(){
  _event=GetEvent();
  if(!IsDATA){
    lhes=GetLHEs();
    gens=GetGens();
    if(IsDYSample){
      GetDYLHEParticles(lhes,lhe_p0,lhe_p1,lhe_l0,lhe_l1,lhe_j0);
      GetDYGenParticles(gens,gen_p0,gen_p1,gen_l0,gen_l1,3);
      GetDYGenParticles(gens,gen_p0,gen_p1,gen_l0_dressed,gen_l1_dressed,1);
      GetDYGenParticles(gens,gen_p0,gen_p1,gen_l0_bare,gen_l1_bare,0);
    }
  }
}
void SMPAnalyzerCore::executeEventWithParameter(Parameter p){
  if(p.weightbit&NominalWeight) FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"lumi",p.w.lumiweight);
  
  if(!_event.PassTrigger(p.triggers)) return;
  if(p.weightbit&NominalWeight) FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"passTrig",p.w.lumiweight);

  if(p.weightbit&NominalWeight){
    FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"PU",p.w.lumiweight*p.w.PUweight);
    FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"prefire",p.w.lumiweight*p.w.PUweight*p.w.prefireweight);
    FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"zpt",p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight);
    FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"z0",p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight);
  }

  double eventweight=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.zptweight;
  if(p.weightbit&NominalWeight) FillHist(p.prefix+p.hprefix+"nlepton"+p.suffix,p.muons.size()+p.electrons.size(),eventweight,10,0,10);

  /////////////////////// selection ///////////////////////
  if(!PassSelection(p)) return;

  ///////////////// efficiency scale factors ///////////////////
  EvalIDSF(p);
  EvalTriggerSF(p);
  if(p.weightbit&NominalWeight){
    FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"RECOSF",eventweight*p.w.electronRECOSF);
    FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"IDSF",eventweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF);
    FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"ISOSF",eventweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF);
    FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"triggerSF",eventweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF);
    FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"CFSF",eventweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF);
  }

  ////// Fill histograms //////////
  FillHists(p);
}
void SMPAnalyzerCore::EvalIDSF(Parameter& p){
  if(!IsDATA){
    if(p.weightbit&EfficiencyWeight){
      p.w.electronRECOSF_sys=fEff->GetStructure(p.k.electronRECOSF);
      p.w.electronIDSF_sys=fEff->GetStructure(p.k.electronIDSF);
      p.w.muonIDSF_sys=fEff->GetStructure(p.k.muonIDSF);
      p.w.muonISOSF_sys=fEff->GetStructure(p.k.muonISOSF);
    }
    for(const auto& electron:p.electrons){
      p.w.electronRECOSF*=fEff->GetEfficiencySF(p.k.electronRECOSF,&electron,0,0);
      if(p.weightbit&EfficiencyWeight){
	int nset=p.w.electronRECOSF_sys.size();
	for(int s=0;s<nset;s++){
	  int nmem=p.w.electronRECOSF_sys[s].size();
	  for(int m=0;m<nmem;m++){
	    p.w.electronRECOSF_sys[s][m]*=fEff->GetEfficiencySF(p.k.electronRECOSF,&electron,s,m);
	  }
	}
      }
      p.w.electronIDSF*=fEff->GetEfficiencySF(p.k.electronIDSF,&electron,0,0);
      if(p.weightbit&EfficiencyWeight){
	int nset=p.w.electronIDSF_sys.size();
	for(int s=0;s<nset;s++){
	  int nmem=p.w.electronIDSF_sys[s].size();
	  for(int m=0;m<nmem;m++){
	    p.w.electronIDSF_sys[s][m]*=fEff->GetEfficiencySF(p.k.electronIDSF,&electron,s,m);
	  }
	}
      }
    }
    for(const auto& muon:p.muons){
      p.w.muonIDSF*=fEff->GetEfficiencySF(p.k.muonIDSF,&muon,0,0);
      if(p.weightbit&EfficiencyWeight){
	int nset=p.w.muonIDSF_sys.size();
	for(int s=0;s<nset;s++){
	  int nmem=p.w.muonIDSF_sys[s].size();
	  for(int m=0;m<nmem;m++){
	    p.w.muonIDSF_sys[s][m]*=fEff->GetEfficiencySF(p.k.muonIDSF,&muon,s,m);
	  }
	}
      }
      p.w.muonISOSF*=fEff->GetEfficiencySF(p.k.muonISOSF,&muon,0,0);
      if(p.weightbit&EfficiencyWeight){
	int nset=p.w.muonISOSF_sys.size();
	for(int s=0;s<nset;s++){
	  int nmem=p.w.muonISOSF_sys[s].size();
	  for(int m=0;m<nmem;m++){
	    p.w.muonISOSF_sys[s][m]*=fEff->GetEfficiencySF(p.k.muonISOSF,&muon,s,m);
	  }
	}
      }
    }
  }
  p.w.CFSF=GetCFSF(p,0);
  p.w.CFSF_up=GetCFSF(p,1);
  p.w.CFSF_down=GetCFSF(p,-1);
}
void SMPAnalyzerCore::EvalTriggerSF(Parameter& p){
  if(!IsDATA){
    vector<Lepton*> triggerables;
    if(p.k.triggerSF.size()){
      if(p.weightbit&EfficiencyWeight){
	p.w.triggerSF_sys=fEff->GetStructure(p.k.triggerSF[0]);
      }
      if(p.k.triggerSF.at(0).Contains("Mu")) triggerables=MakeLeptonPointerVector(p.muons);
      else if(p.k.triggerSF.at(0).Contains("Ele")) triggerables=MakeLeptonPointerVector(p.electrons);
    }
    if(p.k.triggerSF.size()==1){
      p.w.triggerSF*=GetLeptonTriggerSF(p.k.triggerSF[0],triggerables,0,0);
      if(p.weightbit&EfficiencyWeight){
	int nset=p.w.triggerSF_sys.size();
	for(int s=0;s<nset;s++){
	  int nmem=p.w.triggerSF_sys[s].size();
	  for(int m=0;m<nmem;m++){
	    p.w.triggerSF_sys[s][m]*=GetLeptonTriggerSF(p.k.triggerSF[0],triggerables,s,m);
	  }
	}
      }
    }else if(p.k.triggerSF.size()==2){
      if(GetPtThreshold(p.k.triggerSF[0])<GetPtThreshold(p.k.triggerSF[1])){
	p.w.triggerSF*=GetLeptonTriggerORSF(p.k.triggerSF[0],p.k.triggerSF[1],triggerables,0,0);
	if(p.weightbit&EfficiencyWeight){
	  int nset=p.w.triggerSF_sys.size();
	  for(int s=0;s<nset;s++){
	    int nmem=p.w.triggerSF_sys[s].size();
	    for(int m=0;m<nmem;m++){
	      p.w.triggerSF_sys[s][m]*=GetLeptonTriggerORSF(p.k.triggerSF[0],p.k.triggerSF[1],triggerables,s,m);
	    }
	  }
	}
      }else{
	p.w.triggerSF*=GetDileptonTriggerSF(p.k.triggerSF[0],p.k.triggerSF[1],triggerables,0,0);
	if(p.weightbit&EfficiencyWeight){
	  int nset=p.w.triggerSF_sys.size();
	  for(int s=0;s<nset;s++){
	    int nmem=p.w.triggerSF_sys[s].size();
	    for(int m=0;m<nmem;m++){
	      p.w.triggerSF_sys[s][m]*=GetDileptonTriggerSF(p.k.triggerSF[0],p.k.triggerSF[1],triggerables,s,m);
	    }
	  }
	}
      }
    }
  }
}
bool SMPAnalyzerCore::PassSelection(Parameter& p){
  double weight=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.zptweight;

  if(!PassMETFilter()) return false;
  if(p.weightbit&NominalWeight) FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"METfilter",weight);  
  
  if(!p.lepton0||!p.lepton1) return false;
  if(p.weightbit&NominalWeight) FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"Dilepton",weight);

  if(p.c.lepton0pt>0){
    if(p.lepton0->Pt()<p.c.lepton0pt) return false;
  }
  if(p.c.lepton1pt>0){
    if(p.lepton1->Pt()<p.c.lepton1pt) return false;
  }
  if(p.c.muon0pt>0){
    if(p.muons.size()<1) return false;
    if(p.muons.at(0).Pt()<p.c.muon0pt) return false;
  }
  if(p.c.muon1pt>0){
    if(p.muons.size()<2) return false;
    if(p.muons.at(1).Pt()<p.c.muon1pt) return false;
  }
  if(p.c.amuon0pt>0){
    if(p.amuons.size()<1) return false;
    if(p.amuons.at(0).Pt()<p.c.amuon0pt) return false;
  }
  if(p.c.amuon1pt>0){
    if(p.amuons.size()<2) return false;
    if(p.amuons.at(1).Pt()<p.c.amuon1pt) return false;
  }
  if(p.c.electron0pt>0){
    if(p.electrons.size()<1) return false;
    if(p.electrons.at(0).Pt()<p.c.electron0pt) return false;
  }
  if(p.c.electron1pt>0){
    if(p.electrons.size()<2) return false;
    if(p.electrons.at(1).Pt()<p.c.electron1pt) return false;
  }
  if(p.c.aelectron0pt>0){
    if(p.aelectrons.size()<1) return false;
    if(p.aelectrons.at(0).Pt()<p.c.aelectron0pt) return false;
  }
  if(p.c.aelectron1pt>0){
    if(p.aelectrons.size()<2) return false;
    if(p.aelectrons.at(1).Pt()<p.c.aelectron1pt) return false;
  }
  if(p.weightbit&NominalWeight) FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"LepPtCut",weight);

  if(p.lepton0->Charge()*p.lepton1->Charge()>0) p.hprefix+="ss_";
  if(p.weightbit&NominalWeight) FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"charge",weight);

  return true;
}

TH4D* SMPAnalyzerCore::GetHist4D(TString histname){
  TH4D *h = NULL;
  std::map<TString, TH4D*>::iterator mapit = maphist_TH4D.find(histname);
  if(mapit != maphist_TH4D.end()) return mapit->second;
  return h;
}

void SMPAnalyzerCore::FillHist(TString histname,
                            Double_t value_x, Double_t value_y, Double_t value_z, Double_t value_u,
                            Double_t weight,
                            Int_t n_binx, Double_t x_min, Double_t x_max,
                            Int_t n_biny, Double_t y_min, Double_t y_max,
                            Int_t n_binz, Double_t z_min, Double_t z_max,
                            Int_t n_binu, Double_t u_min, Double_t u_max){

  TH4D *this_hist = GetHist4D(histname);
  if( !this_hist ){
    this_hist = new TH4D(histname, "", n_binx, x_min, x_max, n_biny, y_min, y_max, n_binz, z_min, z_max, n_binu, u_min, u_max);
    this_hist->SetDirectory(NULL);
    maphist_TH4D[histname] = this_hist;
  }

  this_hist->Fill(value_x, value_y, value_z, value_u, weight);

}

void SMPAnalyzerCore::FillHist(TString histname,
                            Double_t value_x, Double_t value_y, Double_t value_z, Double_t value_u,
                            Double_t weight,
                            Int_t n_binx, const Double_t *xbins,
                            Int_t n_biny, const Double_t *ybins,
                            Int_t n_binz, const Double_t *zbins,
                            Int_t n_binu, const Double_t *ubins){

  TH4D *this_hist = GetHist4D(histname);
  if( !this_hist ){
    this_hist = new TH4D(histname, "", n_binx, xbins, n_biny, ybins, n_binz, zbins, n_binu, ubins);
    this_hist->SetDirectory(NULL);
    maphist_TH4D[histname] = this_hist;
  }

  this_hist->Fill(value_x, value_y, value_z, value_u, weight);

}

void SMPAnalyzerCore::FillHist(TString histname,
                            Double_t value_x, Double_t value_y, Double_t value_z, Double_t value_u,
                            Double_t weight,
                            Int_t n_binx, const Double_t *xbins,
                            Int_t n_biny, const Double_t *ybins,
                            Int_t n_binz, const Double_t *zbins,
			    Int_t n_binu, Double_t u_min, Double_t u_max){

  TH4D *this_hist = GetHist4D(histname);
  if( !this_hist ){
    TAxis uaxis(n_binu,u_min,u_max);
    vector<double> ubins={};
    for(int i=1;i<n_binu+2;i++) ubins.push_back(uaxis.GetBinLowEdge(i));
    this_hist = new TH4D(histname, "", n_binx, xbins, n_biny, ybins, n_binz, zbins, n_binu, &ubins[0]);
    this_hist->SetDirectory(NULL);
    maphist_TH4D[histname] = this_hist;
  }

  this_hist->Fill(value_x, value_y, value_z, value_u, weight);

}
void SMPAnalyzerCore::FillHists(Parameter& p){
  double weight=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF;
  TLorentzVector dilepton=(*p.lepton0)+(*p.lepton1);
  double dimass=dilepton.M();
  if(dimass>=60&&dimass<120){
    FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"m60to120",weight);
    FillHist(p.prefix+"m60to120/"+p.hprefix+"dimass"+p.suffix,dimass,weight,60,60,120);
    if(dimass>=80&&dimass<100){
      FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"m80to100",weight);
      FillHist(p.prefix+"m80to100/"+p.hprefix+"dimass"+p.suffix,dimass,weight,40,80,100);
    }
  }
}

void SMPAnalyzerCore::WriteHist(){
  AnalyzerCore::WriteHist();
  outfile->cd();
  for(std::map< TString, TH4D* >::iterator mapit = maphist_TH4D.begin(); mapit!=maphist_TH4D.end(); mapit++){
    TString this_fullname=mapit->second->GetName();
    TString this_name=this_fullname(this_fullname.Last('/')+1,this_fullname.Length());
    TString this_suffix=this_fullname(0,this_fullname.Last('/'));
    TDirectory *dir = outfile->GetDirectory(this_suffix);
    if(!dir){
      outfile->mkdir(this_suffix);
    }
    outfile->cd(this_suffix);
    mapit->second->Write(this_name);
    outfile->cd();
  }
}

void SMPAnalyzerCore::FillHist(TString histname, double value, map<TString,double> weights, int n_bin, double x_min, double x_max){
  for(const auto& [suffix,weight]:weights) FillHist(histname+suffix,value,weight,n_bin,x_min,x_max);
}
void SMPAnalyzerCore::FillHist(TString histname, double value, map<TString,double> weights, int n_bin, double *xbins){
  for(const auto& [suffix,weight]:weights) FillHist(histname+suffix,value,weight,n_bin,xbins);
}
void SMPAnalyzerCore::FillHist(TString histname, double value_x, double value_y, map<TString,double> weights, int n_binx, double x_min, double x_max, int n_biny, double y_min, double y_max){
  for(const auto& [suffix,weight]:weights) FillHist(histname+suffix,value_x,value_y,weight,n_binx,x_min,x_max,n_biny,y_min,y_max);
}
void SMPAnalyzerCore::FillHist(TString histname, double value_x, double value_y, map<TString,double> weights, int n_binx, double *xbins, int n_biny, double *ybins){
  for(const auto& [suffix,weight]:weights) FillHist(histname+suffix,value_x,value_y,weight,n_binx,xbins,n_biny,ybins);
}
void SMPAnalyzerCore::FillHist(TString histname, double value_x, double value_y, double value_z, map<TString,double> weights, int n_binx, double x_min, double x_max, int n_biny, double y_min, double y_max, int n_binz, double z_min, double z_max){
  for(const auto& [suffix,weight]:weights) FillHist(histname+suffix,value_x,value_y,value_z,weight,n_binx,x_min,x_max,n_biny,y_min,y_max,n_binz,z_min,z_max);
}
void SMPAnalyzerCore::FillHist(TString histname, double value_x, double value_y, double value_z, map<TString,double> weights, int n_binx, double *xbins, int n_biny, double *ybins, int n_binz, double *zbins){
  for(const auto& [suffix,weight]:weights) FillHist(histname+suffix,value_x,value_y,value_z,weight,n_binx,xbins,n_biny,ybins,n_binz,zbins);
}
void SMPAnalyzerCore::FillHist(TString histname, double value_x, double value_y, double value_z, double value_u, map<TString,double> weights, int n_binx, double x_min, double x_max, int n_biny, double y_min, double y_max, int n_binz, double z_min, double z_max, int n_binu, double u_min, double u_max){
  for(const auto& [suffix,weight]:weights) FillHist(histname+suffix,value_x,value_y,value_z,value_u,weight,n_binx,x_min,x_max,n_biny,y_min,y_max,n_binz,z_min,z_max,n_binu,u_min,u_max);
}
void SMPAnalyzerCore::FillHist(TString histname, double value_x, double value_y, double value_z, double value_u, map<TString,double> weights, int n_binx, double *xbins, int n_biny, double *ybins, int n_binz, double *zbins, int n_binu, double *ubins){
  for(const auto& [suffix,weight]:weights) FillHist(histname+suffix,value_x,value_y,value_z,value_u,weight,n_binx,xbins,n_biny,ybins,n_binz,zbins,n_binu,ubins);
}
void SMPAnalyzerCore::FillHist(TString histname, double value_x, double value_y, double value_z, double value_u, map<TString,double> weights, int n_binx, double *xbins, int n_biny, double *ybins, int n_binz, double *zbins, int n_binu, double u_min, double u_max){
  for(const auto& [suffix,weight]:weights) FillHist(histname+suffix,value_x,value_y,value_z,value_u,weight,n_binx,xbins,n_biny,ybins,n_binz,zbins,n_binu,u_min,u_max);
}

void SMPAnalyzerCore::FillDileptonHists(TString pre,TString suf,Particle *l0,Particle *l1,double w){
  TLorentzVector dilepton=*l0+*l1;
  double dimass=dilepton.M();
  double dipt=dilepton.Pt();
  double dirap=dilepton.Rapidity();
  FillHist(pre+"dimass"+suf,dimass,w,200,0,400);
  FillHist(pre+"dipt"+suf,dipt,w,200,0,400);
  FillHist(pre+"dirap"+suf,dirap,w,120,-6,6);
  vector<Particle*> leps;
  if(l0->Pt()>l1->Pt()) leps={l0,l1};
  else leps={l1,l0};
  for(int i=0;i<(int)leps.size();i++){
    FillHist(Form("%sl%dpt%s",pre.Data(),i,suf.Data()),leps.at(i)->Pt(),w,200,0,400);
    FillHist(Form("%sl%deta%s",pre.Data(),i,suf.Data()),leps.at(i)->Eta(),w,100,-5,5);
    FillHist(Form("%sl%dphi%s",pre.Data(),i,suf.Data()),leps.at(i)->Phi(),w,80,-4,4);
  }
  FillHist(pre+"lldelR"+suf,l0->DeltaR(*l1),w,70,0,7);  
  FillHist(pre+"lldelphi"+suf,l0->DeltaPhi(*l1),w,80,-4,4);
  FillHist(pre+"lldeleta"+suf,fabs(l0->Eta()-l1->Eta()),w,100,-5,5);
}
double SMPAnalyzerCore::GetPtThreshold(TString path){
  TString str=path(TRegexp("[0-9]+"));
  if(str.IsFloat()) return str.Atof();
  else return -1.;
}
bool SMPAnalyzerCore::IsExists(TString filepath){
  ifstream fcheck(filepath);
  return fcheck.good();
}
void SMPAnalyzerCore::SetupEfficiency(){
  TString configpath=getenv("DATA_DIR")+TString("/")+GetEra()+"/ID/eff.conf";
  if(IsExists(configpath)){
    fEff=new EfficiencyTool(configpath);
  }
}
void SMPAnalyzerCore::DeleteEfficiency(){
  if(fEff) delete fEff;
}
double SMPAnalyzerCore::GetLeptonTriggerSF(TString triggerSF_key,const vector<Lepton*>& leps,int set,int mem){
  if(IsDATA) return 1;
  if(triggerSF_key=="") return 1;
  if(triggerSF_key=="Default") return 1;

  double data_eff=1.,sim_eff=1.;
  for(const auto& lep:leps){
    data_eff*=1-fEff->GetDataEfficiency(triggerSF_key,lep,set,mem);
    sim_eff*=1-fEff->GetSimEfficiency(triggerSF_key,lep,set,mem);
  }
  data_eff=1-data_eff;
  sim_eff=1-sim_eff;
  if(sim_eff==0) return 1.;
  else return data_eff/sim_eff;
}
double SMPAnalyzerCore::GetLeptonTriggerORSF(TString triggerSF_key0,TString triggerSF_key1,const vector<Lepton*>& leps,int set,int mem){
  if(IsDATA) return 1;

  double lumi0=1.; //only trigger0 on
  double lumi1=0.; //only trigger1 on
  double lumi2=0.; //both off
  if(DataYear==2017&&triggerSF_key0.Contains("IsoMu24")&&triggerSF_key1.Contains("IsoMu27")){
    lumi0=37997.005;    lumi1=3480.873;    lumi2=0.;
  }else if(DataYear==2017&&triggerSF_key0.Contains("Ele27")&&triggerSF_key1.Contains("Ele32")){
    lumi0=31661.026;    lumi1=9522.208;    lumi2=0.295;
  }else if(DataYear==2018&&triggerSF_key0.Contains("Ele28")&&triggerSF_key1.Contains("Ele32")){
    lumi0=23687.253;    lumi1=36140.626;   lumi2=0.;
  }else{
    cout<<"[SMPAnalyzerCore::LeptonTriggerOR_SF] not available combination "<<triggerSF_key0<<"||"<<triggerSF_key1<<" for "<<DataEra<<endl;
    exit(EXIT_FAILURE);
  }

  double data_eff_key0=1.,data_eff_key1=1.,sim_eff=1.;
  for(const auto& lep:leps){
    data_eff_key0*=1-fEff->GetDataEfficiency(triggerSF_key0,lep,set,mem);
    data_eff_key1*=1-fEff->GetDataEfficiency(triggerSF_key1,lep,set,mem);
    sim_eff*=1-fEff->GetSimEfficiency(triggerSF_key0,lep,set,mem);
  }
  data_eff_key0=1-data_eff_key0;
  data_eff_key1=1-data_eff_key1;
  sim_eff=1-sim_eff;
  if(sim_eff==0) return 1.;
  else return (lumi0*data_eff_key0+lumi1*data_eff_key1)/((lumi0+lumi1+lumi2)*sim_eff);
}
double SMPAnalyzerCore::GetDileptonTriggerSF(TString triggerSF_key0,TString triggerSF_key1,const vector<Lepton*>& leps,int set,int mem){
  if(IsDATA) return 1;
  if((triggerSF_key0==""||triggerSF_key0=="Default")&&(triggerSF_key1==""||triggerSF_key1=="Default")) return 1;
  int nlep=leps.size();
  if(nlep<2){
    cout<<"[SMPAnalyzerCore::DileptonTrigger_SF] nlep < 2. return 1."<<endl;
    return 1.;
  }
  double data_noleg1=1.,sim_noleg1=1.;
  vector<double> data_oneleg1_noleg2(nlep,1.);
  vector<double> sim_oneleg1_noleg2(nlep,1.);
  for(int i=0;i<nlep;i++){
    double data_eff_leg1=fEff->GetDataEfficiency(triggerSF_key0,leps.at(i),set,mem);
    double data_eff_leg2=fEff->GetDataEfficiency(triggerSF_key1,leps.at(i),set,mem);
    double sim_eff_leg1=fEff->GetSimEfficiency(triggerSF_key0,leps.at(i),set,mem);
    double sim_eff_leg2=fEff->GetSimEfficiency(triggerSF_key1,leps.at(i),set,mem);
    data_noleg1*=(1-data_eff_leg1);
    sim_noleg1*=(1-sim_eff_leg1);
    for(int j=0;j<nlep;j++){
      if(i==j){
	data_oneleg1_noleg2[j]*=data_eff_leg1;
	sim_oneleg1_noleg2[j]*=sim_eff_leg1;
      }else{
	data_oneleg1_noleg2[j]*=(1-data_eff_leg2);
	sim_oneleg1_noleg2[j]*=(1-sim_eff_leg2);
      }
    }
  }
  double data_eff=1.-data_noleg1;
  double sim_eff=1.-sim_noleg1;
  for(int i=0;i<nlep;i++){
    data_eff-=data_oneleg1_noleg2[i];
    sim_eff-=sim_oneleg1_noleg2[i];
  }
  if(sim_eff==0) return 1.;
  else return data_eff/sim_eff;
}
void SMPAnalyzerCore::SetupZptWeight(){
  /// TODO: currently, temp Zpt weight from 2017 double muon for all other channels
  cout<<"[SMPAnalyzerCore::SetupZptWeight] setting zptcor"<<endl;
  TString datapath=getenv("DATA_DIR");
  if(!IsExists(datapath+"/"+GetEra()+"/SMP/ZptWeight.root")){
    cout<<"[SMPAnalyzerCore::SetupZptWeight] no ZptWeight.root"<<endl;
    return;
  }
  TFile fzpt(datapath+"/"+GetEra()+"/SMP/ZptWeight.root");
  for(const auto&& key:*(fzpt.GetListOfKeys())){
    TH2D* this_hist=(TH2D*)((TKey*)key)->ReadObj();
    TString histname=this_hist->GetName();
    TString striphistname=histname;
    int extent=0,start=striphistname.Index(TRegexp("_iter[0-9]*$"),&extent);
    if(start>=0) striphistname.Replace(start,extent,"");
    cout<<"[SMPAnalyzerCore::SetupZptWeight] setting "<<histname<<endl;
    if(histname.Contains(TRegexp("_iter0$"))){
      map_hist_zpt[striphistname]=this_hist;
      this_hist->SetDirectory(0);
    }else if(histname.Contains(TRegexp("_iter[0-9]*$"))){
      map_hist_zpt[striphistname]->Multiply(this_hist);
    }else if(histname.Contains(TRegexp("_norm$"))){
      map_hist_zpt[histname]=this_hist;
      this_hist->SetDirectory(0);
    }
  }
}
void SMPAnalyzerCore::SetupRoccoR(){
  cout<<"[SMPAnalyzerCore::SetupRoccoR] setting Rocheseter Correction"<<endl;
  TString datapath=getenv("DATA_DIR");
  TString rocpath=datapath+"/"+GetEra()+"/RoccoR/RoccoR"+GetEraShort()+"UL.txt";
  if(IsExists(rocpath)) roc=new RoccoR(rocpath.Data());
  else cout<<"[SMPAnalyzerCore::SetupRoccoR] no "+rocpath<<endl;
  TString erashort=GetEraShort();
  TString rocelepath=datapath+"/"+GetEra()+"/RoccoR/e_"+erashort(2,3)+"UL.txt";
  if(DataYear==2016) rocelepath=datapath+"/"+GetEra()+"/RoccoR/e_"+erashort(2,3)+"UL_1.txt";
  if(IsExists(rocelepath)) rocele=new Aepcor(rocelepath.Data());
  else cout<<"[SMPAnalyzerCore::SetupRoccoR] no "+rocelepath<<endl;  
}
void SMPAnalyzerCore::SetupZ0Weight(){
  cout<<"[SMPAnalyzerCore::SetupZ0Weight] setting Z0Weight"<<endl;
  TString datapath=getenv("DATA_DIR");
  if(!IsExists(datapath+"/"+GetEra()+"/SMP/Z0Weight.root")){
    cout<<"[SMPAnalyzerCore::SetupZ0Weight] no Z0Weight.root"<<endl;
    return;
  }
  TFile fz0(datapath+"/"+GetEra()+"/SMP/Z0Weight.root");
  hz0_data=(TF1*)fz0.Get("data_fit");
  hz0_mc=(TF1*)fz0.Get("mc_fit");
  fz0.Close();
}
double SMPAnalyzerCore::GetZ0Weight(double valx){
  return 1.; /// FIXME: no Z0 weight for UL
  double data_val = hz0_data->Eval(valx);
  double mc_val = hz0_mc->Eval(valx);
  double norm = hz0_mc->Integral(-100,100)/hz0_data->Integral(-100,100);
  return norm*data_val/mc_val;
}
double SMPAnalyzerCore::GetZptWeight(double zpt,double zrap,Lepton::Flavour flavour){
  double valzptcor=1.;
  double valzptcor_norm=1.;
  zrap=fabs(zrap);
  TString sflavour=flavour==Lepton::MUON?"muon":"electron";
  TString MCName=MCSample;
  if(MCName.Contains(TRegexp("^DY[0-9]Jets$"))) MCName="DYJets";
  if(MCName=="DYJets_Summer20") MCName="DYJets";
  if(MCName.Contains(TRegexp("^DYJets_Pt-[0-9]*To[0-9Inf]*$"))) MCName="DYJets";
  if(MCName.Contains(TRegexp("^DYJets_M-[0-9]*to[0-9Inf]*$"))) MCName="DYJets";
  TString hzptname=MCName+"_"+sflavour;
  TString hnormname=hzptname+"_norm";
  auto it=map_hist_zpt.find(hzptname);
  if(it!=map_hist_zpt.end()){
    valzptcor*=GetBinContentUser(map_hist_zpt[hzptname],zpt,zrap,0);
  }
  it=map_hist_zpt.find(hnormname);
  if(it!=map_hist_zpt.end()){
    valzptcor_norm*=GetBinContentUser(map_hist_zpt[hnormname],zpt,zrap,0);
  }
  return valzptcor*valzptcor_norm;
}
void SMPAnalyzerCore::SetupCFRate(){
  cout<<"[SMPAnalyzerCore::SetupCFRate] setting CFRate"<<endl;
  TString datapath=getenv("DATA_DIR");
  if(!IsExists(datapath+"/"+GetEra()+"/SMP/CFRate.root")){
    cout<<"[SMPAnalyzerCore::SetupZ0Weight] no CFRate.root"<<endl;
    return;
  }
  TFile f(datapath+"/"+GetEra()+"/SMP/CFRate.root");
  hcfrate_data=(TH2*)f.Get("cfdata");
  if(hcfrate_data){
    cout<<"[SMPAnalyzerCore::SetupCFRate] load hcfrate_data"<<endl;
    hcfrate_data->SetDirectory(0);
  }else cout<<"[SMPAnalyzerCore::SetupCFRate] no hcfrate_data"<<endl;
  hcfrate_mc=(TH2*)f.Get("cfmc");
  if(hcfrate_mc){
    cout<<"[SMPAnalyzerCore::SetupCFRate] load hcfrate_mc"<<endl;
    hcfrate_mc->SetDirectory(0);
  }else cout<<"[SMPAnalyzerCore::SetupCFRate] no hcfrate_mc"<<endl;
  hcfsf=(TH2*)f.Get("cfsf");
  if(hcfsf){
    cout<<"[SMPAnalyzerCore::SetupCFRate] load hcfsf"<<endl;
    hcfsf->SetDirectory(0);
  }else cout<<"[SMPAnalyzerCore::SetupCFRate] no hcfsf"<<endl;
  f.Close();
}
double SMPAnalyzerCore::GetCFSF(const Lepton* l,int sys){
  if(IsDATA) return 1.;
  if(!hcfsf) return 1.;
  if(!l) return 1.;
  if(l->LeptonFlavour()!=Lepton::ELECTRON) return 1.;
  return GetBinContentUser(hcfsf,l->Eta(),l->Pt(),sys);
}
double SMPAnalyzerCore::GetCFSF(const Parameter& p,int sys){
  if(IsDATA) return 1.;
  double sf=1.;
  if(p.lepton0&&!p.truth_lepton0.IsEmpty())
    if(p.lepton0->Charge()*p.truth_lepton0.Charge()<0) sf*=GetCFSF(p.lepton0,sys);
  if(p.lepton1&&!p.truth_lepton1.IsEmpty())
    if(p.lepton1->Charge()*p.truth_lepton1.Charge()<0) sf*=GetCFSF(p.lepton1,sys);
  return sf;
}
void SMPAnalyzerCore::DeleteCFRate(){
  if(hcfrate_data) delete hcfrate_data;
  if(hcfrate_mc) delete hcfrate_mc;
  if(hcfsf) delete hcfsf;
}
void SMPAnalyzerCore::PrintGens(const vector<Gen>& gens){
  cout<<"index\tpid\tmother\tstatus\tpropt\thard\n";
  for(int i=0;i<(int)gens.size();i++){
    gens[i].Print();
    cout<<gens.at(i).Index()<<"\t"<<gens.at(i).PID()<<"\t"<<gens.at(i).MotherIndex()<<"\t"<<gens.at(i).Status()<<"\t"<<gens.at(i).isPrompt()<<"\t"<<gens.at(i).isHardProcess()<<endl;
  }
}

double SMPAnalyzerCore::GetBinContentUser(TH2* hist,double valx,double valy,int sys){
  double xmin=hist->GetXaxis()->GetXmin();
  double xmax=hist->GetXaxis()->GetXmax();
  double ymin=hist->GetYaxis()->GetXmin();
  double ymax=hist->GetYaxis()->GetXmax();
  if(xmin>=0) valx=fabs(valx);
  if(valx<xmin) valx=xmin+0.001;
  if(valx>=xmax) valx=xmax-0.001;
  if(ymin>=0) valy=fabs(valy);
  if(valy<ymin) valy=ymin+0.001;
  if(valy>=ymax) valy=ymax-0.001;
  return hist->GetBinContent(hist->FindBin(valx,valy))+sys*hist->GetBinError(hist->FindBin(valx,valy));
}

double SMPAnalyzerCore::GetBinContentUser(TH3* hist,double valx,double valy,double valz,int sys){
  double xmin=hist->GetXaxis()->GetXmin();
  double xmax=hist->GetXaxis()->GetXmax();
  double ymin=hist->GetYaxis()->GetXmin();
  double ymax=hist->GetYaxis()->GetXmax();
  double zmin=hist->GetZaxis()->GetXmin();
  double zmax=hist->GetZaxis()->GetXmax();
  if(xmin>=0) valx=fabs(valx);
  if(valx<xmin) valx=xmin+0.001;
  if(valx>xmax) valx=xmax-0.001;
  if(ymin>=0) valy=fabs(valy);
  if(valy<ymin) valy=ymin+0.001;
  if(valy>ymax) valy=ymax-0.001;
  if(zmin>=0) valz=fabs(valz);
  if(valz<zmin) valz=zmin+0.001;
  if(valz>zmax) valz=zmax-0.001;
  return hist->GetBinContent(hist->FindBin(valx,valy,valz))+sys*hist->GetBinError(hist->FindBin(valx,valy,valz));
}
void SMPAnalyzerCore::GetDYLHEParticles(const vector<LHE>& lhes,LHE& p0,LHE& p1,LHE& l0,LHE& l1,LHE& j0){
  if(!IsDYSample){
    cout <<"[AFBAnalyzer::GetDYLHEParticles] this is for DY event"<<endl;
    exit(EXIT_FAILURE);
  }
  p0=LHE();
  p1=LHE();
  l0=LHE();
  l1=LHE();
  j0=LHE();
  for(int i=0;i<(int)lhes.size();i++){
    if(p0.ID()==0&&lhes[i].Status()==-1&&lhes[i].Eta()>0) p0=lhes[i];
    if(p1.ID()==0&&lhes[i].Status()==-1&&lhes[i].Eta()<0) p1=lhes[i];
    if(l0.ID()==0&&(abs(lhes[i].ID())==11||abs(lhes[i].ID())==13||abs(lhes[i].ID())==15)) l0=lhes[i];
    if(l0.ID()&&lhes[i].ID()==-l0.ID()) l1=lhes[i];
    if(lhes[i].Status()==1)
      if(abs(lhes[i].ID())<=6||lhes[i].ID()==21)
	if(lhes[i].Pt()>j0.Pt()) j0=lhes[i];
  }
  if(p0.ID()==0||p1.ID()==0||l0.ID()==0||l1.ID()==0){
    cout <<"[AFBAnalyzer::GetLHEParticles] something is wrong"<<endl;
    exit(EXIT_FAILURE);
  }
  if(l0.Pt()<l1.Pt()){
    LHE temp=l0;
    l0=l1;
    l1=temp;
  }
}


void SMPAnalyzerCore::GetDYGenParticles(const vector<Gen>& gens,Gen& parton0,Gen& parton1,Gen& l0,Gen& l1,int mode){
  //mode 0:bare 1:dressed01 2:dressed04 3:beforeFSR
  if(!IsDYSample){
    cout <<"[SMPAnalyzerCore::GetDYGenParticles] this is for DY event"<<endl;
    exit(EXIT_FAILURE);
  }
  parton0=Gen();
  parton1=Gen();
  l0=Gen();
  l1=Gen();
  vector<const Gen*> leptons;
  vector<const Gen*> photons;
  int ngen=gens.size();
  for(int i=0;i<ngen;i++){
    if(!gens.at(i).isPrompt()) continue;
    int genpid=gens.at(i).PID();
    if(gens.at(i).isHardProcess()){
      if(abs(genpid)<7||genpid==21){
	if(parton0.IsEmpty()) parton0=gens[i];
	else if(parton1.IsEmpty()) parton1=gens[i];
      }
    }
    if(gens.at(i).Status()==1){
      if(abs(genpid)==11||abs(genpid)==13) leptons.push_back(&gens[i]);
      else if(gens.at(i).PID()==22) photons.push_back(&gens[i]);
    }
  }
  int nlepton=leptons.size();
  for(int i=0;i<nlepton;i++){
    for(int j=i+1;j<nlepton;j++){
      if(!(leptons[i]->PID()+leptons[j]->PID()==0)) continue;
      if((*leptons[i]+*leptons[j]).M()>(l0+l1).M()){
	if(leptons[i]->Pt()>leptons[j]->Pt()){
	  l0=*leptons[i];
	  l1=*leptons[j];
	}else{
	  l0=*leptons[j];
	  l1=*leptons[i];
	}
      }
    }
  }
  if(mode>=3){
    if(nlepton>=4){
      for(int i=0;i<nlepton;i++){
	if(leptons[i]->Index()==l0.Index()||leptons[i]->Index()==l1.Index()) continue;
	for(int j=i+1;j<nlepton;j++){
	  if(leptons[j]->Index()==l0.Index()||leptons[j]->Index()==l1.Index()) continue;
	  if(!(leptons[i]->PID()+leptons[j]->PID()==0)) continue;
	  vector<int> history_i=TrackGenSelfHistory(*leptons[i],gens);
	  vector<int> history_j=TrackGenSelfHistory(*leptons[j],gens);
	  if(history_i.at(1)==history_j.at(1)){
	    photons.push_back(leptons[i]);
	    photons.push_back(leptons[j]);
	  }
	}
      }
    }	
    for(const auto& photon:photons){
      vector<int> history=TrackGenSelfHistory(*photon,gens);
      if(gens[history.at(1)].PID()==l0.PID()) l0+=*photon;
      else if(gens[history.at(1)].PID()==l1.PID()) l1+=*photon;
      else if(gens[history.at(1)].PID()==23){ // for minnlo+photos
	if(photon->DeltaR(l0)<photon->DeltaR(l1)) l0+=*photon;
	else l1+=*photon;
      }
    }    
  }else if(mode>=1){
    double delr=mode==1?0.1:0.4;
    for(const auto& photon:photons){
      if(l0.DeltaR(*photon)>delr&&l1.DeltaR(*photon)>delr) continue;
      if(l0.DeltaR(*photon)<l1.DeltaR(*photon)) l0+=*photon;
      else l1+=*photon;
    }
  }
}

Gen SMPAnalyzerCore::SMPGetGenMatchedLepton(const Lepton& lep,const std::vector<Gen>& gens,int mode){
  //0: default
  //1: dressed 0.1
  Gen gen_lepton=GetGenMatchedLepton(lep,gens);
  if(gen_lepton.IsEmpty()) return gen_lepton;
  if(mode==1){ // dressed 0.1 cone
    for(const auto& gen: gens){
      if(gen.Status()!=1) continue;
      if(gen.PID()!=22) continue;
      if(gen.DeltaR(gen_lepton)>0.1) continue;
      gen_lepton+=gen;
    }
  }
    
  return gen_lepton;
}

std::vector<Electron> SMPAnalyzerCore::SMPGetElectrons(TString id, double ptmin, double fetamax){
  std::vector<Electron> out;
  if(id=="passMediumID_Selective"){
    std::vector<Electron> electrons = GetAllElectrons();
    for(unsigned int i=0; i<electrons.size(); i++){
      Electron this_electron= electrons.at(i);
      if(!( this_electron.Pt()>ptmin ))	continue;
      if(!( fabs(this_electron.scEta())<fetamax )) continue;
      if(!( this_electron.PassID("passMediumID") ))	continue;
      if(!electron_isGsfCtfScPixChargeConsistent->at(i)) continue;
      out.push_back(this_electron);
    }
  }else if(id=="passTightID_Selective"){
    std::vector<Electron> electrons = GetAllElectrons();
    for(unsigned int i=0; i<electrons.size(); i++){
      Electron this_electron= electrons.at(i);
      if(!( this_electron.Pt()>ptmin ))	continue;
      if(!( fabs(this_electron.scEta())<fetamax )) continue;
      if(!( this_electron.PassID("passTightID") )) continue;
      if(!electron_isGsfCtfScPixChargeConsistent->at(i)) continue;
      out.push_back(this_electron);
    }
  }else if(id=="passMediumIDWithAntiIso"){
    vector<Electron> electrons = GetAllElectrons();
    for(unsigned int i=0; i<electrons.size(); i++){
      Electron el= electrons.at(i);
      if(!( el.Pt()>ptmin ))	continue;
      if(!( fabs(el.scEta())<fetamax )) continue;
      if( fabs(el.scEta()) <= 1.479 ){
	if(! (el.Full5x5_sigmaIetaIeta() < 0.0106) ) continue;
	if(! (fabs(el.dEtaSeed()) < 0.0032) ) continue;
	if(! (fabs(el.dPhiIn()) < 0.0547) ) continue;
	if( (el.HoverE() < 0.046 + 1.16/el.scE() + 0.0324*el.Rho()/el.scE()) && (el.RelIso() < 0.0478+0.506/el.UncorrPt()) ) continue;
	if(! (fabs(el.InvEminusInvP()) < 0.184) ) continue;
	if(! (el.NMissingHits() <= 1) ) continue;
	if(! (el.PassConversionVeto()) ) continue;
      }else{
	if(! (el.Full5x5_sigmaIetaIeta() < 0.0387) ) continue;
	if(! (fabs(el.dEtaSeed()) < 0.00632) ) continue;
	if(! (fabs(el.dPhiIn()) <  0.0394 ) ) continue;
	if( (el.HoverE() < 0.0275 + 2.52/el.scE() + 0.183*el.Rho()/el.scE()) && (el.RelIso() < 0.0658+0.963/el.UncorrPt()) ) continue;
	if(! (fabs(el.InvEminusInvP()) < 0.0721) ) continue;
	if(! (el.NMissingHits() <= 1) ) continue;
	if(! (el.PassConversionVeto()) ) continue;
      }
      out.push_back(el);
    }
  }else if(id=="passAntiMediumID"){
    vector<Electron> electrons = GetAllElectrons();
    for(unsigned int i=0; i<electrons.size(); i++){
      Electron el= electrons.at(i);
      if(!( el.Pt()>ptmin ))	continue;
      if(!( fabs(el.scEta())<fetamax )) continue;
      if( el.PassID("passMediumID") ) continue;
      out.push_back(el);
    }
  }else if(id=="passAntiLooseID"){
    vector<Electron> electrons = GetAllElectrons();
    for(unsigned int i=0; i<electrons.size(); i++){
      Electron el= electrons.at(i);
      if(!( el.Pt()>ptmin ))	continue;
      if(!( fabs(el.scEta())<fetamax )) continue;
      if( el.PassID("passLooseID") ) continue;
      out.push_back(el);
    }
  }else out=GetElectrons(id,ptmin,fetamax);
  std::sort(out.begin(),out.end(),PtComparing);
  return out;
}    
std::vector<Muon> SMPAnalyzerCore::SMPGetMuons(TString id,double ptmin,double fetamax){
  vector<Muon> out;
  if(id=="POGTightWithLooseTrkIso"){
    vector<Muon> muons=GetMuons("POGTight",ptmin,fetamax);
    for(auto const& muon: muons){
      //if(muon.TrkIso()/muon.Pt()<0.1) out.push_back(muon);
      if(muon.PassSelector(Muon::Selector::TkIsoLoose)) out.push_back(muon);
    }
  }else if(id=="POGMediumWithLooseTrkIso"){
    vector<Muon> muons=GetMuons("POGMedium",ptmin,fetamax);
    for(auto const& muon: muons){
      if(muon.PassSelector(Muon::Selector::TkIsoLoose)) out.push_back(muon);
    }
  }else if(id=="POGMediumWithAntiLooseTrkIso"){
    vector<Muon> muons=GetMuons("POGMedium",ptmin,fetamax);
    for(auto const& muon: muons){
      if(!muon.PassSelector(Muon::Selector::TkIsoLoose)) out.push_back(muon);
    }
  }else if(id=="POGTightWithAntiIso"){
    vector<Muon> muons=GetMuons("POGTight",ptmin,fetamax);
    for(auto const& muon: muons){
      if(muon.RelIso()>0.3) out.push_back(muon);
    }
  }else if(id=="POGTightWithAntiMediumIso"){
    vector<Muon> muons=GetMuons("POGTight",ptmin,fetamax);
    for(auto const& muon: muons){
      if(muon.RelIso()>0.2) out.push_back(muon);
    }
  }else out=GetMuons(id,ptmin,fetamax);
  return out;
}

std::vector<Muon> SMPAnalyzerCore::MuonMomentumCorrection(const vector<Muon>& muons,int sys,int set,int member){
  if(!roc) return std::vector<Muon>(muons);
  std::vector<Muon> out;
  for(auto muon:muons){
    double rc=1.;
    double rcerr=0.;
    if(set>=0){
      if(IsDATA){
	rc=roc->kScaleDT(muon.Charge(),muon.MiniAODPt(),muon.Eta(),muon.Phi(),set,member);
	rcerr=roc->kScaleDTerror(muon.Charge(),muon.MiniAODPt(),muon.Eta(),muon.Phi());
      }else{
	Gen gen=GetGenMatchedLepton(muon,gens);
	if(gen.IsEmpty()){
	  double u=gRandom->Rndm();
	  rc=roc->kSmearMC(muon.Charge(),muon.MiniAODPt(),muon.Eta(),muon.Phi(),muon.TrackerLayers(),u,set,member);
	  rcerr=roc->kSmearMCerror(muon.Charge(),muon.MiniAODPt(),muon.Eta(),muon.Phi(),muon.TrackerLayers(),u);
	}else{
	  rc=roc->kSpreadMC(muon.Charge(),muon.MiniAODPt(),muon.Eta(),muon.Phi(),gen.Pt(),set,member);
	  rcerr=roc->kSpreadMCerror(muon.Charge(),muon.MiniAODPt(),muon.Eta(),muon.Phi(),gen.Pt());
	}
      }      
    }
    muon.SetPtEtaPhiM(muon.MiniAODPt()*(rc+sys*rcerr),muon.Eta(),muon.Phi(),muon.M());
    out.push_back(muon);
  }
  std::sort(out.begin(),out.end(),PtComparing);
  return out;
}

std::vector<Electron> SMPAnalyzerCore::ElectronEnergyCorrection(const vector<Electron>& electrons,int set,int member){
  if(!rocele) return std::vector<Electron>(electrons);
  std::vector<Electron> out;
  double u=gRandom->Rndm();
  for(auto electron:electrons){
    if(set>=0){
      double rc=1.;
      //double rcerr=0.;
      double el_eta=electron.scEta();
      double el_phi=electron.Phi();
      // FIXME: rocelecor don't have the corrections for 2.4<|eta|<2.5
      if(el_eta>=2.4) el_eta=2.39;
      if(el_eta<=-2.4) el_eta=-2.39;
      if(fabs(el_eta)<2.4){
	if(IsDATA){
	  rc=rocele->kScaleDT(electron.UncorrPt(),el_eta,el_phi,electron.R9(),run,set,member);
	  electron*=rc*electron.UncorrE()/electron.E();
	}else{	
	  Gen gen=SMPGetGenMatchedLepton(electron,gens,1);
	  if(!gen.IsEmpty()&&gen.Pt()/electron.Pt()>0.5){
	    rc=rocele->kSpreadMC(electron.UncorrPt(),el_eta,el_phi,electron.R9(),u,gen.Pt(),set,member);
	  }else{
	    rc=rocele->kSmearMC(electron.UncorrPt(),el_eta,el_phi,electron.R9(),u,set,member);
	  }
	  if(TMath::IsNaN(rc)) rc=1.;
	  electron*=rc*electron.UncorrE()/electron.E();
	}      
      }
    }else if(set==-1){ //no energe cor
      electron*=electron.UncorrE()/electron.E();
    }else{
      cout<<"[SMPAnalyzerCore::ElectronEnergyCorrection] wrong set "<<set<<endl;
      exit(ENODATA);
    }
    out.push_back(electron);
  }
  std::sort(out.begin(),out.end(),PtComparing);
  return out;
}
  
void SMPAnalyzerCore::FillCutflow(TString histname,TString label,double weight){
  TH1D* hist=NULL;
  auto it=maphist_TH1D.find(histname);
  if(it==maphist_TH1D.end()){
    hist=new TH1D(histname,"",1,0,1);
    hist->GetXaxis()->SetBinLabel(1,label);
    maphist_TH1D[histname]=hist;    
  }else hist=it->second;

  int nbin=hist->GetNbinsX();
  int ibin=0;
  for(int i=1;i<=nbin;i++){
    if(hist->GetXaxis()->GetBinLabel(i)==label){
      ibin=i;
    }
  }

  if(!ibin){
    hist->SetBins(nbin+1,0,nbin+1);
    ibin=nbin+1;
    hist->GetXaxis()->SetBinLabel(ibin,label);
  }
  hist->Fill(ibin-0.5,weight);
}
    
TString SMPAnalyzerCore::Replace(TString str,TRegexp reg,TString repl){
  int extent;
  int start=str.Index(reg,&extent);
  if(start>=0) return str.Replace(start,extent,repl);
  else return str;
}

SMPAnalyzerCore::Parameter::Parameter(){
}
SMPAnalyzerCore::Parameter::~Parameter(){
}
void SMPAnalyzerCore::Parameter::SetChannel(TString ch){
  vector<TString> availables={"el","ee","eE","EE","mu","mm","mM","MM","em","me"};
  bool pass=false;
  for(const TString& avail:availables)
    if(ch==avail) pass=true;
  if(!pass){
    cout<<"[SMPAnalyzerCore::Parameter::SetChannel] not available channel "<<ch<<endl;
    exit(EXIT_FAILURE);
  }
  channel=ch;
  SetLeptons();
}
void SMPAnalyzerCore::Parameter::SetElectronKeys(TString elID,vector<TString> trig){
  k.electronRECOSF="Electron_RECO";
  k.electronIDSF=elID;
  k.triggerSF=trig;
}
void SMPAnalyzerCore::Parameter::SetMuonKeys(TString muID,TString muISO,vector<TString> trig){
  k.muonIDSF=muID;
  k.muonISOSF=muISO;
  k.triggerSF=trig;
}
void SMPAnalyzerCore::Parameter::SetLeptonPtCut(double l0pt,double l1pt){
  c.lepton0pt=l0pt;
  c.lepton1pt=l1pt;
}
void SMPAnalyzerCore::Parameter::SetLeptons(){
  leptons={};
  lepton0=NULL; 
  lepton1=NULL;
  truth_lepton0=Gen();
  truth_lepton1=Gen();
  if(channel=="") return;
  unsigned int ie=0,iae=0,im=0,iam=0;
  int nc=channel.Length();
  for(int ic=0;ic<nc;ic++){
    char c=channel[ic];
    if(c=='e'||c=='l'){
      if(electrons.size()>ie){
	leptons.push_back(&electrons.at(ie));
	ie++;
      }else leptons.push_back(NULL);
    }else if(c=='m'||c=='u'){
      if(muons.size()>im){
	leptons.push_back(&muons.at(im));
	im++;
      }else leptons.push_back(NULL);
    }else if(c=='E'){
      if(aelectrons.size()>iae){
	leptons.push_back(&aelectrons.at(iae));
	iae++;
      }else leptons.push_back(NULL);
    }else if(c=='M'){
      if(amuons.size()>iam){
	leptons.push_back(&amuons.at(iam));
	iam++;
      }else leptons.push_back(NULL);
    }
  }
  if(leptons.size()>0) lepton0=leptons.at(0);
  if(leptons.size()>1) lepton1=leptons.at(1);
  if(lepton0) truth_lepton0=SMPGetGenMatchedLepton(*lepton0,gens);
  if(lepton1) truth_lepton1=SMPGetGenMatchedLepton(*lepton1,gens);
}
void SMPAnalyzerCore::Parameter::SetGens(vector<Gen> gs){
  gens=gs;
  SetLeptons();
}
void SMPAnalyzerCore::Parameter::SetElectrons(vector<Electron> els){
  electrons=els;
  SetLeptons();
}    
void SMPAnalyzerCore::Parameter::SetMuons(vector<Muon> mus){
  muons=mus;
  SetLeptons();
}    
void SMPAnalyzerCore::Parameter::SetAElectrons(vector<Electron> els){
  aelectrons=els;
  SetLeptons();
}    
void SMPAnalyzerCore::Parameter::SetAMuons(vector<Muon> mus){
  amuons=mus;
  SetLeptons();
}    

SMPAnalyzerCore::Parameter SMPAnalyzerCore::MakeParameter(TString channel,TString option){
  Parameter p;
  p.SetChannel(channel);
  p.option=option;
  p.SetGens(gens);

  p.hprefix="";
  p.w.lumiweight=reductionweight;
  p.w.PUweight=1;  p.w.PUweight_up=1;  p.w.PUweight_down=1;
  p.w.prefireweight=1;  p.w.prefireweight_up=1;  p.w.prefireweight_down=1;
  p.w.z0weight=1;
  p.w.zptweight=1;
  if(!IsDATA){
    p.w.lumiweight*=weight_norm_1invpb*_event.MCweight()*_event.GetTriggerLumi("Full");
    p.w.PUweight=mcCorr->GetPileUpWeight(nPileUp,0);
    p.w.PUweight_up=mcCorr->GetPileUpWeight(nPileUp,1);
    p.w.PUweight_down=mcCorr->GetPileUpWeight(nPileUp,-1);
    p.w.z0weight=GetZ0Weight(vertex_Z);
    if(DataYear<2018){
      p.w.prefireweight=L1PrefireReweight_Central;
      p.w.prefireweight_up=L1PrefireReweight_Up;
      p.w.prefireweight_down=L1PrefireReweight_Down;
    }
    if(IsDYSample){
      if(abs(lhe_l0.ID())==11||abs(lhe_l0.ID())==13){
	TLorentzVector genZ=(gen_l0+gen_l1);
	p.w.zptweight=GetZptWeight(genZ.Pt(),genZ.Rapidity(),abs(lhe_l0.ID())==13?Lepton::Flavour::MUON:Lepton::Flavour::ELECTRON);
      }else p.hprefix+="tau_";
    }
  }

  p.prefix=p.channel+GetEraShort()+"/";
  if(p.channel=="mu"){
    p.SetMuonKeys("Muon_MediumID_trkIsoLoose","",{"IsoMu24_MediumID_trkIsoLoose"});
    p.SetMuons(MuonMomentumCorrection(SMPGetMuons("POGMediumWithLooseTrkIso",0.0,2.4),0,3));
    p.SetLeptonPtCut(27,10);
    if(GetEraShort()=="2016a"){
      p.triggers={"HLT_IsoMu24_v","HLT_IsoTkMu24_v"};
    }else if(GetEraShort()=="2016b"){
      p.triggers={"HLT_IsoMu24_v","HLT_IsoTkMu24_v"};
    }else if(GetEraShort()=="2017"){
      p.triggers={"HLT_IsoMu24_v","HLT_IsoMu27_v"};
      p.k.triggerSF={"IsoMu24_MediumID_trkIsoLoose","IsoMu27_MediumID_trkIsoLoose"};
    }else if(GetEraShort()=="2018"){
      p.triggers={"HLT_IsoMu24_v"};
    }
  }else if(p.channel=="mm"){
    p.SetMuonKeys("Muon_MediumID_trkIsoLoose","",{"Mu17Leg1_MediumID_trkIsoLoose","Mu8Leg2_MediumID_trkIsoLoose"});
    p.SetMuons(MuonMomentumCorrection(SMPGetMuons("POGMediumWithLooseTrkIso",0.0,2.4),0,3));
    p.SetLeptonPtCut(20,10);
    if(GetEraShort()=="2016a"){
      p.triggers={"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v","HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v","HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v","HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",};
    }else if(GetEraShort()=="2016b"){
      p.triggers={"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v","HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v","HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v","HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",};
    }else if(GetEraShort()=="2017"){
      p.triggers={"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v"};
    }else if(GetEraShort()=="2018"){
      p.triggers={"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v"};
    }
  }else if(p.option.Contains("TightIDSelectiveCharge")&&p.channel=="el"){
    p.prefix="tight/"+p.prefix;
    p.SetElectronKeys("Electron_TightID_Selective",{"Ele27_TightID_Selective"});
    p.SetElectrons(ElectronEnergyCorrection(SMPGetElectrons("passTightID_Selective",0.0,2.4),0,0));
    //p.SetElectrons(SMPGetElectrons("passTightID_Selective",0.0,2.4));
    p.SetLeptonPtCut(30,10);
    if(GetEraShort()=="2016a"){
      p.triggers={"HLT_Ele27_WPTight_Gsf_v"};
    }else if(GetEraShort()=="2016b"){
      p.triggers={"HLT_Ele27_WPTight_Gsf_v"};
    }else if(GetEraShort()=="2017"){
      p.triggers={"HLT_Ele27_WPTight_Gsf_v","HLT_Ele32_WPTight_Gsf_v"};
      p.k.triggerSF={"Ele27_TightID_Selective","Ele32_TightID_Selective"};
    }else if(GetEraShort()=="2018"){
      p.triggers={"HLT_Ele28_WPTight_Gsf_v","HLT_Ele32_WPTight_Gsf_v"};
      p.k.triggerSF={"Ele28_TightID_Selective","Ele32_TightID_Selective"};
    }
  }else if(p.channel=="el"){
    p.SetElectronKeys("Electron_MediumID",{"Ele27_MediumID"});
    p.SetElectrons(ElectronEnergyCorrection(SMPGetElectrons("passMediumID",0.0,2.4),0,0));
    //p.SetElectrons(SMPGetElectrons("passMediumID",0.0,2.4));
    p.SetLeptonPtCut(30,10);
    if(GetEraShort()=="2016a"){
      p.triggers={"HLT_Ele27_WPTight_Gsf_v"};
    }else if(GetEraShort()=="2016b"){
      p.triggers={"HLT_Ele27_WPTight_Gsf_v"};
    }else if(GetEraShort()=="2017"){
      p.triggers={"HLT_Ele27_WPTight_Gsf_v","HLT_Ele32_WPTight_Gsf_v"};
      p.k.triggerSF={"Ele27_MediumID","Ele32_MediumID"};
    }else if(GetEraShort()=="2018"){
      p.triggers={"HLT_Ele28_WPTight_Gsf_v","HLT_Ele32_WPTight_Gsf_v"};
      p.k.triggerSF={"Ele28_MediumID","Ele32_MediumID"};
    }
  }else if(p.channel=="ee"){
    p.SetElectronKeys("Electron_MediumID",{"Ele23Leg1_MediumID","Ele12Leg2_MediumID"});
    p.SetElectrons(ElectronEnergyCorrection(SMPGetElectrons("passMediumID",0.0,2.4),0,0));
    //p.SetElectrons(SMPGetElectrons("passMediumID",0.0,2.4));
    p.SetLeptonPtCut(25,15);
    if(GetEraShort()=="2016a") p.triggers={"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"};
    else if(GetEraShort()=="2016b") p.triggers={"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"};
    else if(GetEraShort()=="2017") p.triggers={"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v"};
    else if(GetEraShort()=="2018") p.triggers={"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v"};
  }else if(p.channel=="me"){
    p.SetMuonKeys("Muon_MediumID_trkIsoLoose","",{"IsoMu24_MediumID_trkIsoLoose"});
    p.k.electronIDSF="Electron_MediumID";
    p.SetMuons(MuonMomentumCorrection(SMPGetMuons("POGMediumWithLooseTrkIso",0.0,2.4),0,3));
    p.SetElectrons(SMPGetElectrons("passMediumID",0.0,2.4));
    p.SetLeptonPtCut(27,10);
    if(GetEraShort()=="2016a"){
      p.triggers={"HLT_IsoMu24_v","HLT_IsoTkMu24_v"};
    }else if(GetEraShort()=="2016b"){
      p.triggers={"HLT_IsoMu24_v","HLT_IsoTkMu24_v"};
    }else if(GetEraShort()=="2017"){
      p.triggers={"HLT_IsoMu24_v","HLT_IsoMu27_v"};
      p.k.triggerSF={"IsoMu24_MediumID_trkIsoLoose","IsoMu27_MediumID_trkIsoLoose"};
    }else if(GetEraShort()=="2018"){
      p.triggers={"HLT_IsoMu24_v"};
    }
  }else if(p.channel=="em"){
    p.SetElectronKeys("Electron_MediumID",{"Ele27_MediumID"});
    p.k.muonIDSF="Muon_MediumID_trkIsoLoose";
    p.SetElectrons(SMPGetElectrons("passMediumID",0.0,2.4));
    p.SetMuons(MuonMomentumCorrection(SMPGetMuons("POGMediumWithLooseTrkIso",0.0,2.4),0,3));
    p.SetLeptonPtCut(30,10);
    if(GetEraShort()=="2016a"){
      p.triggers={"HLT_Ele27_WPTight_Gsf_v"};
    }else if(GetEraShort()=="2016b"){
      p.triggers={"HLT_Ele27_WPTight_Gsf_v"};
    }else if(GetEraShort()=="2017"){
      p.triggers={"HLT_Ele27_WPTight_Gsf_v","HLT_Ele32_WPTight_Gsf_v"};
      p.k.triggerSF={"Ele27_MediumID","Ele32_MediumID"};
    }else if(GetEraShort()=="2018"){
      p.triggers={"HLT_Ele28_WPTight_Gsf_v","HLT_Ele32_WPTight_Gsf_v"};
      p.k.triggerSF={"Ele28_MediumID","Ele32_MediumID"};
    }
  }else if(p.channel=="mM"){
    p.SetMuonKeys("Muon_MediumID_trkIsoLoose","",{"IsoMu24_MediumID_trkIsoLoose"});
    p.SetMuons(MuonMomentumCorrection(SMPGetMuons("POGMediumWithLooseTrkIso",0.0,2.4),0,3));
    p.SetAMuons(MuonMomentumCorrection(SMPGetMuons("POGMediumWithAntiLooseTrkIso",0.0,2.4),0,3));
    p.SetLeptonPtCut(27,10);
    if(GetEraShort()=="2016a"){
      p.triggers={"HLT_IsoMu24_v","HLT_IsoTkMu24_v"};
    }else if(GetEraShort()=="2016b"){
      p.triggers={"HLT_IsoMu24_v","HLT_IsoTkMu24_v"};
    }else if(GetEraShort()=="2017"){
      p.triggers={"HLT_IsoMu24_v","HLT_IsoMu27_v"};
      p.k.triggerSF={"IsoMu24_MediumID_trkIsoLoose","IsoMu27_MediumID_trkIsoLoose"};
    }else if(GetEraShort()=="2018"){
      p.triggers={"HLT_IsoMu24_v"};
    }
  }else if(p.channel=="MM"){
    p.k.muonIDSF="Muon_MediumID_trkIsoLoose";
    p.SetAMuons(MuonMomentumCorrection(SMPGetMuons("POGMediumWithAntiLooseTrkIso",0.0,2.4),0,3));
    p.SetLeptonPtCut(20,10);
    if(GetEraShort()=="2016a"){
      p.triggers={"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v","HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v","HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v","HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",};
    }else if(GetEraShort()=="2016b"){
      p.triggers={"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v","HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v","HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v","HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",};
    }else if(GetEraShort()=="2017"){
      p.triggers={"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v"};
    }else if(GetEraShort()=="2018"){
      p.triggers={"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v"};
    }
  }else if(p.channel=="eE"){
    p.SetElectronKeys("Electron_MediumID",{"Ele27_MediumID"});
    p.SetElectrons(SMPGetElectrons("passMediumID",0.0,2.4));
    p.SetAElectrons(SMPGetElectrons("passAntiLooseID",0.0,2.4));
    p.SetLeptonPtCut(30,10);
    if(GetEraShort()=="2016a"){
      p.triggers={"HLT_Ele27_WPTight_Gsf_v"};
    }else if(GetEraShort()=="2016b"){
      p.triggers={"HLT_Ele27_WPTight_Gsf_v"};
    }else if(GetEraShort()=="2017"){
      p.triggers={"HLT_Ele27_WPTight_Gsf_v","HLT_Ele32_WPTight_Gsf_v"};
      p.k.triggerSF={"Ele27_MediumID","Ele32_MediumID"};
    }else if(GetEraShort()=="2018"){
      p.triggers={"HLT_Ele28_WPTight_Gsf_v","HLT_Ele32_WPTight_Gsf_v"};
      p.k.triggerSF={"Ele28_MediumID","Ele32_MediumID"};
    }
  }else if(p.channel=="EE"){
    p.k.electronIDSF="Electron_MediumID";
    p.SetAElectrons(SMPGetElectrons("passAntiLooseID",0.0,2.4));
    p.SetLeptonPtCut(25,15);
    if(GetEraShort()=="2016a") p.triggers={"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"};
    else if(GetEraShort()=="2016b") p.triggers={"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"};
    else if(GetEraShort()=="2017") p.triggers={"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v"};
    else if(GetEraShort()=="2018") p.triggers={"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v"};
  }
  return p;
}
vector<TString> SMPAnalyzerCore::Split(TString s,TString del){
  TObjArray* array=s.Tokenize(del);
  vector<TString> out;
  for(const auto& obj:*array){
    out.push_back(((TObjString*)obj)->String());
  }
  array->Delete();
  return out;
}
