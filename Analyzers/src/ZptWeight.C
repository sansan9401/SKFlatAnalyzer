#include "ZptWeight.h"
const double ZptWeight::massbin[ZptWeight::massbinnum+1]={52,60,70,80,90,100,120,150,200,300,400,1000};
const double ZptWeight::ybin[ZptWeight::ybinnum+1]={0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0};
const double ZptWeight::ptbin[ZptWeight::ptbinnum+1]={0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,3.25,3.5,3.75,4,4.25,4.5,4.75,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18,18.5,19,19.5,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,42,44,46,48,50,52,54,56,59,62,66,70,75,80,90,100,120,140,160,180,200,240,280,320,360,400,450,500,550,650};
ZptWeight::ZptWeight(){
}
ZptWeight::~ZptWeight(){
}
void ZptWeight::executeEvent(){
  ////////////////////////check genlevel//////////////////
  if(IsDYSample){
    if(abs(lhe_l0.ID())!=15){
      Parameter p;
      if(abs(lhe_l0.ID())==11) p=MakeParameter("ee");
      else if(abs(lhe_l1.ID())==13) p=MakeParameter("mm");
      TLorentzVector genZ=(gen_l0+gen_l1);
      FillHist(p.prefix+"gen_mypt_nozptweight",genZ.M(),fabs(genZ.Rapidity()),genZ.Pt(),p.w.lumiweight,massbinnum,massbin,ybinnum,ybin,ptbinnum,ptbin);
      FillHist(p.prefix+"gen_mypt",genZ.M(),fabs(genZ.Rapidity()),genZ.Pt(),p.w.lumiweight*p.w.zptweight,massbinnum,massbin,ybinnum,ybin,ptbinnum,ptbin);
    }
  }

  if(!IsDATA||DataStream.Contains("SingleMuon")) 
    executeEventWithParameter(MakeParameter("mu"));
  if(!IsDATA||DataStream.Contains("DoubleMuon")) 
    executeEventWithParameter(MakeParameter("mm"));
  if(!IsDATA||DataStream.Contains("SingleElectron")||DataStream.Contains("EGamma")) 
    executeEventWithParameter(MakeParameter("el"));
  if(!IsDATA||DataStream.Contains("DoubleEG")||DataStream.Contains("EGamma")) 
    executeEventWithParameter(MakeParameter("ee"));
}

void ZptWeight::FillHists(Parameter& p){
  TLorentzVector dilepton=(*p.lepton0)+(*p.lepton1);
  double dimass=dilepton.M();
  double dipt=dilepton.Pt();
  double dirap=fabs(dilepton.Rapidity());
  TLorentzVector genZ=(gen_l0+gen_l1);

  map<TString,double> weightmap;
  weightmap[""]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF;
  weightmap["_nozptweight"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF;

  for(const auto& [wname,w]:weightmap){
    TString pre=p.prefix+p.hprefix;
    TString suf=p.suffix+wname;    
    FillHist(pre+"myptgpt"+suf,dimass,dirap,dipt,genZ.Pt(),w,massbinnum,massbin,ybinnum,ybin,ptbinnum,ptbin,ptbinnum,ptbin);
  }
}
