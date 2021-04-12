#include "ZptWeight.h"

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
      if(genZ.M()>80&&genZ.M()<100){
	FillHist(p.prefix+"gen_diptdirap",genZ.Pt(),fabs(genZ.Rapidity()),p.w.lumiweight*p.w.zptweight,ptbinnum,ptbin,genrapbinnum,genrapbin);
	FillHist(p.prefix+"gen_diptdirap_nozptweight",genZ.Pt(),fabs(genZ.Rapidity()),p.w.lumiweight,ptbinnum,ptbin,genrapbinnum,genrapbin);
	FillHist(p.prefix+"gen_dipt2dirap",pow(genZ.Pt(),2),fabs(genZ.Rapidity()),p.w.lumiweight*p.w.zptweight,pt2binnum,pt2bin,genrapbinnum,genrapbin);
	FillHist(p.prefix+"gen_dipt2dirap_nozptweight",pow(genZ.Pt(),2),fabs(genZ.Rapidity()),p.w.lumiweight,pt2binnum,pt2bin,genrapbinnum,genrapbin);
      }
      FillHist(p.prefix+"gen_diptdimass",genZ.Pt(),genZ.M(),p.w.lumiweight*p.w.zptweight,ptbinnum,ptbin,massbinnum,massbin);
      FillHist(p.prefix+"gen_diptdimass_nozptweight",genZ.Pt(),genZ.M(),p.w.lumiweight,ptbinnum,ptbin,massbinnum,massbin);
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
  double dirap=dilepton.Rapidity();
  if(dimass<80||dimass>100) return;

  map<TString,double> weightmap;
  weightmap[""]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF;
  weightmap["_nozptweight"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF;

  for(const auto& [wname,w]:weightmap){
    TString pre=p.prefix+"m80to100/"+p.hprefix;
    TString suf=p.suffix+wname;
    
    FillHist(pre+"diptdirap"+suf,dipt,fabs(dirap),w,ptbinnum,ptbin,rapbinnum,rapbin);
    FillHist(pre+"dipt"+suf,dipt,w,100,0,200);
    for(int i=0;i<rapbinnum;i++){
      if(fabs(dirap)>rapbin[i]&&fabs(dirap)<rapbin[i+1]){
	FillHist(pre+Form("dipt_bin%d",i)+suf,dipt,w,100,0,200);
	break;
      }
    }
  }
}
