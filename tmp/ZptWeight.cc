TString analyzer="ZptWeight";
TString recoskim="SkimTree_Dilepton";
TString genskim="";

TString dyname="DYJetsToMuMu_MiNNLO";
TString taudyname="DYJetsToTauTau_MiNNLO";
vector<TString> backgrounds={"TTLL_powheg","WJets_MG","WW_pythia","WZ_pythia","ZZ_pythia"};
map<TString,vector<TString>> datas={
  {"ee2016", {"DoubleEG_B_ver2","DoubleEG_C","DoubleEG_D","DoubleEG_E","DoubleEG_F","DoubleEG_G","DoubleEG_H"}},
  {"mm2016", {"DoubleMuon_B_ver2","DoubleMuon_C","DoubleMuon_D","DoubleMuon_E","DoubleMuon_F","DoubleMuon_G","DoubleMuon_H"}},
  {"2016",{"DoubleEG_B_ver2","DoubleEG_C","DoubleEG_D","DoubleEG_E","DoubleEG_F","DoubleEG_G","DoubleEG_H","DoubleMuon_B_ver2","DoubleMuon_C","DoubleMuon_D","DoubleMuon_E","DoubleMuon_F","DoubleMuon_G","DoubleMuon_H"}},
  {"mm2016b", {"DoubleMuon_F","DoubleMuon_G","DoubleMuon_H"}},
  {"ee2017", {"DoubleEG_B","DoubleEG_C","DoubleEG_D","DoubleEG_E","DoubleEG_F"}},
  {"mm2017", {"DoubleMuon_B","DoubleMuon_C","DoubleMuon_D","DoubleMuon_E","DoubleMuon_F"}},
  {"mu2017", {"SingleMuon_B","SingleMuon_C","SingleMuon_D","SingleMuon_E","SingleMuon_F"}},
  {"2017",{"DoubleEG_B","DoubleEG_C","DoubleEG_D","DoubleEG_E","DoubleEG_F","DoubleMuon_B","DoubleMuon_C","DoubleMuon_D","DoubleMuon_E","DoubleMuon_F"}},
  {"ee2018", {"EGamma_A","EGamma_B","EGamma_C","EGamma_D"}},
  {"mm2018", {"DoubleMuon_A","DoubleMuon_B","DoubleMuon_C","DoubleMuon_D"}},
  {"2018", {"EGamma_A","EGamma_B","EGamma_C","EGamma_D","DoubleMuon_A","DoubleMuon_B","DoubleMuon_C","DoubleMuon_D"}},
};

TH1* GetHist(TString filename,TString histname){
  TFile f(filename);
  TH1* hist=(TH1*)f.Get(histname);
  if(!hist){
    cout<<"No "<<histname<<" in "<<filename<<endl;
    return hist;
  }
  hist->SetDirectory(0);
  return hist;
}


TH1* GetZptWeight(TString mode){
  TH1::SetDefaultSumw2(true);
  TString prefix=mode;
  TString syear=mode(2,4);
  if(mode.Contains("2016a")) syear="2016preVFP";
  else if(mode.Contains("2016b")) syear="2016postVFP";
  TString schannel;
  if(mode.BeginsWith("mm")) schannel="muon";
  else if(mode.BeginsWith("mu")) schannel="muon";
  else if(mode.BeginsWith("ee")) schannel="electron";
  else {
    cout<<"Unknown mode "<<mode<<endl;
    exit(1);
  }      
  TString mcfileprefix=TString(getenv("SKFlatOutputDir"))+getenv("SKFlatV")+"/"+analyzer+"/"+syear+"/"+analyzer+"_";
  TString datafileprefix=TString(getenv("SKFlatOutputDir"))+getenv("SKFlatV")+"/"+analyzer+"/"+syear+"/DATA/"+analyzer+"_";
  if(recoskim!=""){
    mcfileprefix+=recoskim+"_";
    datafileprefix+=recoskim+"_";
  }
  
  TH1* hdy=GetHist(mcfileprefix+dyname+".root",prefix+"/m80to100/diptdirap");
  TH1* hdata=NULL;
  for(const auto& data:datas[mode]){
    TString file=datafileprefix+data+".root";
    cout<<"hdata+="<<file<<endl;
    if(hdata){
      hdata->Add(GetHist(file,prefix+"/m80to100/diptdirap"));
    }else{
      hdata=GetHist(file,prefix+"/m80to100/diptdirap");
    }
  }
  cout<<"hdata-= tau from "<<mcfileprefix+dyname+".root"<<endl;
  if(taudyname!="") hdata->Add(GetHist(mcfileprefix+taudyname+".root",prefix+"/m80to100/tau_diptdirap"),-1);
  else hdata->Add(GetHist(mcfileprefix+dyname+".root",prefix+"/m80to100/tau_diptdirap"),-1);
  for(const auto& bkname:backgrounds){
    TString file=mcfileprefix+bkname+".root";
    cout<<"hdata-="<<file<<endl;
    hdata->Add(GetHist(file,prefix+"/m80to100/diptdirap"),-1);
  }
  hdata->Scale(1./hdata->Integral());
  hdy->Scale(1./hdy->Integral());
  hdata->Divide(hdy);
  hdata->SetName(dyname+"_"+schannel);
  return hdata;
}
TH1* GetNormWeight(TString mode){
  TH1::SetDefaultSumw2(true);
  TString schannel;
  TString syear=mode(2,4);
  if(mode.Contains("2016a")) syear="2016preVFP";
  else if(mode.Contains("2016b")) syear="2016postVFP";
  TString prefix=mode;
  if(mode.BeginsWith("mm")) schannel="muon";
  else if(mode.BeginsWith("mu")){
    schannel="muon";
    prefix=mode.ReplaceAll("mu","mm");
  }else if(mode.BeginsWith("ee")) schannel="electron";
  else {
    cout<<"Unknown mode "<<mode<<endl;
    exit(1);
  }
  TString mcfileprefix=TString(getenv("SKFlatOutputDir"))+getenv("SKFlatV")+"/"+analyzer+"/"+syear+"/"+analyzer+"_";
  if(genskim!=""){
    mcfileprefix+=genskim+"_";
  }
  
  TH2D* hdy=(TH2D*)GetHist(mcfileprefix+dyname+".root",prefix+"/gen_diptdirap");
  TH2D* hdy_nozptcor=(TH2D*)GetHist(mcfileprefix+dyname+".root",prefix+"/gen_diptdirap_nozptweight");
  TH1D* dirap=hdy->ProjectionY("dirap");
  TH1D* dirap_nozptcor=hdy_nozptcor->ProjectionY("dirap_nozptweight");
  vector<double> xbin={hdy->GetXaxis()->GetBinLowEdge(1),hdy->GetXaxis()->GetBinLowEdge(hdy->GetNbinsX()+1)};
  vector<double> ybin;
  for(int i=1;i<dirap->GetNbinsX()+2;i++){
    ybin.push_back(dirap->GetXaxis()->GetBinLowEdge(i));
  }
  ybin.push_back(5);
  TH2D* norm=new TH2D(dyname+"_"+schannel+"_norm","",1,&xbin[0],ybin.size()-1,&ybin[0]);
  for(int i=1;i<dirap->GetNbinsX()+2;i++){
    norm->SetBinContent(1,i,dirap_nozptcor->GetBinContent(i)/dirap->GetBinContent(i));
  }
  
  return (TH1*)norm;
}
void SaveZptWeight(TString mode){
  TString syear=mode(2,4);
  if(mode.Contains("2016a")) syear="2016preVFP";
  else if(mode.Contains("2016b")) syear="2016postVFP";
  TString schannel;
  if(mode.BeginsWith("mm")) schannel="muon";
  else if(mode.BeginsWith("mu")) schannel="muon";
  else if(mode.BeginsWith("ee")) schannel="electron";
  else {
    cout<<"Unknown mode "<<mode<<endl;
    exit(1);
  }      
  TString zptfiledir=getenv("SKFlat_WD")+TString("/data/")+getenv("SKFlatV")+"/"+syear+"/SMP/";
  system("mkdir -p "+zptfiledir);
  TFile f(zptfiledir+"ZptWeight.root","update");
  int i=0;
  while(f.Get(dyname+"_"+schannel+Form("_iter%d",i))){
    i++;
  }
  TH1* hist=GetZptWeight(mode);
  f.cd();
  hist->SetName(dyname+"_"+schannel+Form("_iter%d",i));
  hist->Write();
}
void SaveNormWeight(TString mode){
  TString syear=mode(2,4);
  if(mode.Contains("2016a")) syear="2016preVFP";
  else if(mode.Contains("2016b")) syear="2016postVFP";
  TString schannel;
  if(mode.BeginsWith("mm")) schannel="muon";
  else if(mode.BeginsWith("mu")) schannel="muon";
  else if(mode.BeginsWith("ee")) schannel="electron";
  else {
    cout<<"Unknown mode "<<mode<<endl;
    exit(1);
  }      
  TH1* hist=GetNormWeight(mode);
  TString zptfiledir=getenv("SKFlat_WD")+TString("/data/")+getenv("SKFlatV")+"/"+syear+"/SMP/";
  system("mkdir -p "+zptfiledir);
  TFile f(zptfiledir+"ZptWeight.root","update");
  hist->SetName(dyname+"_"+schannel+"_norm");
  hist->Write();
}
    
void Iterate(TString mode,int n=3){
  if(n<1){
    cout<<"n should be greater than 0 "<<endl;
    exit(1);
  }					 
  TString prefix=mode;
  TString syear=mode(2,5);
  TString schannel;
  if(mode.BeginsWith("mm")) schannel="muon";
  else if(mode.BeginsWith("mu")) schannel="muon";
  else if(mode.BeginsWith("ee")) schannel="electron";
  else if(mode.Contains(TRegexp("^201[6-8][ab]?$"))){
    syear=mode;
  }else{
    cout<<"Unknown mode "<<mode<<endl;
    exit(1);
  }      
  TString cmd;
  for(const auto& data:datas[mode]){
    cmd+="SKFlat.py -a ZptWeight -n 10 -e "+syear+" -i "+data(0,data.Index('_'))+" -p "+data(data.Index('_')+1,999);
    if(recoskim!="") cmd+=" --skim "+recoskim;
    cmd+=" &";
  }
  cmd+="SKFlat.py -a ZptWeight -n 40 -e "+syear+" -i "+dyname;
  if(recoskim!="") cmd+=" --skim "+recoskim;
  cmd+=" &";
  if(taudyname!=""){
    cmd+="SKFlat.py -a ZptWeight -n 10 -e "+syear+" -i "+taudyname;
    if(recoskim!="") cmd+=" --skim "+recoskim;
    cmd+=" &";
  }
  for(const auto& bk:backgrounds){
    cmd+="SKFlat.py -a ZptWeight -n 10 -e "+syear+" -i "+bk;
    if(recoskim!="") cmd+=" --skim "+recoskim;
    cmd+=" &";
  }
  cmd+=" wait;";
  cout<<"Iteration 0, "<<mode<<endl;
  cout<<cmd<<endl;
  //system(cmd); //FIXME uncomment this
  cout<<"Save ZptWeight iter0"<<endl;
  if(mode.Contains(TRegexp("^[0-9]*$"))){
    SaveZptWeight("ee"+syear);
    SaveZptWeight("mm"+syear);
  }else{
    SaveZptWeight(mode);
  }    
  for(int i=1;i<n;i++){
    TString dycmd="SKFlat.py -a ZptWeight -n 100 -e "+syear+" -i "+dyname;
    if(recoskim!="") dycmd+=" --skim "+recoskim;
    cout<<"Iteration "<<i<<", "<<mode<<endl;
    cout<<dycmd<<endl;
    system(dycmd);
    cout<<"Save ZptWeight iter"<<i<<endl;
    if(mode.Contains(TRegexp("^[0-9]*$"))){
      SaveZptWeight("ee"+syear);
      SaveZptWeight("mm"+syear);
    }else{
      SaveZptWeight(mode);
    }    
  }
  TString gendycmd="SKFlat.py -a ZptWeight -n 100 -e "+syear+" -i "+dyname;
  if(genskim!="") gendycmd+=" --skim "+genskim;
  cout<<"Iteration for normalization, "<<mode<<endl;
  system(gendycmd);
  cout<<"Save ZptWeight for normalization"<<endl;
  if(mode.Contains(TRegexp("^[0-9]*$"))){
    SaveNormWeight("ee"+syear);
    SaveNormWeight("mm"+syear);
  }else{
    SaveNormWeight(mode);
  }    
  return;
}
