#include "EfficiencyTool.h"

Efficiency::Efficiency(){
}
Efficiency::~Efficiency(){
  vector<vector<vector<TH2*>>*> targets={&fDataPlus,&fDataMinus,&fSimPlus,&fSimMinus};
  for(auto vv:targets){
    for(auto v:*vv){
      for(auto h:v){
	if(h) delete h;
      }
    }
    vv->clear();
  }
}
void Efficiency::AddSet(bool isData,vector<TString> paths,int charge){
  if(charge==0){
    AddSet(isData,paths,+1);
    AddSet(isData,paths,-1);
  }else{
    vector<vector<TH2*>> *target=GetTarget(isData,charge);
    vector<TH2*> hists;
    for(const TString& path:paths){
      TString file=path(0,path.Index(':'));
      TString key=path(path.Index(':')+1,path.Length());
      if(!HasKey(file,key)){
	cout<<"[Efficiency::AddSet] no "<<key<<" in "<<file<<endl;
	exit(ENODATA);
      }
      TFile f(file);
      TH2* hist=(TH2*)f.Get(key);
      hist->SetDirectory(NULL);
      hists.push_back(hist);
    }
    target->push_back(hists);
  }
}
void Efficiency::AddDataSet(vector<TString> paths,int charge){
  AddSet(true,paths,charge);
}
void Efficiency::AddSimSet(vector<TString> paths,int charge){
  AddSet(false,paths,charge);
}
void Efficiency::AddSetReplica(bool isData,TString nominal,TString stat,int nreplica,int charge){
  if(charge==0){
    AddSetReplica(isData,nominal,stat,nreplica,+1);
    AddSetReplica(isData,nominal,stat,nreplica,-1);
  }else{
    vector<vector<TH2*>> *target=GetTarget(isData,charge);
    vector<TH2*> hists;
    TH2* central=NULL;
    {
      TString path=nominal;
      TString file=path(0,path.Index(':'));
      TString key=path(path.Index(':')+1,path.Length());
      if(!HasKey(file,key)){
	cout<<"[Efficiency::AddSetReplica] no "<<key<<" in "<<file<<endl;
	exit(ENODATA);
      }
      TFile f(file);
      central=(TH2*)f.Get(key);
      central->SetDirectory(NULL);
      hists.push_back(central);
    }
    {
      TString path=stat;
      TString file=path(0,path.Index(':'));
      TString key=path(path.Index(':')+1,path.Length());
      if(!HasKey(file,key)){
	cout<<"[Efficiency::AddSetReplica] no "<<key<<" in "<<file<<endl;
	exit(ENODATA);
      }
      TFile f(file);
      TH2* hstat=(TH2*)f.Get(key);
      TRandom3 random(path.Hash());
      if(SameContents(central,hstat)){
	for(int i=1;i<nreplica;i++){
	  TH2* hist=(TH2*)central->Clone();
	  hist->SetDirectory(NULL);
	  for(int b=0;b<hist->GetNcells();b++){
	    double err=hstat->GetBinError(b);
	    double val=random.Gaus(central->GetBinContent(b),err);
	    hist->SetBinContent(b,val);
	    hist->SetBinError(b,0);
	  }
	  hists.push_back(hist);
	}
      }else{ 
	for(int i=1;i<nreplica;i++){
	  TH2* hist=(TH2*)central->Clone();
	  hist->SetDirectory(NULL);
	  for(int b=0;b<hist->GetNcells();b++){
	    double err=hstat->GetBinContent(b);
	    double val=random.Gaus(central->GetBinContent(b),err);
	    if(val>1) val=1;
	    if(val<0) val=0;
	    hist->SetBinContent(b,val);
	    hist->SetBinError(b,0);
	  }
	  hists.push_back(hist);
	}      
      }
    }
    target->push_back(hists);
  }
}
void Efficiency::AddDataSetReplica(TString nominal,TString stat,int nreplica,int charge){
  AddSetReplica(true,nominal,stat,nreplica,charge);
}
void Efficiency::AddSimSetReplica(TString nominal,TString stat,int nreplica,int charge){
  AddSetReplica(false,nominal,stat,nreplica,charge);
}
void Efficiency::AddSetUpDown(bool isData,TString path,int charge){
  if(charge==0){
    AddSetUpDown(isData,path,+1);
    AddSetUpDown(isData,path,-1);
  }else{
    vector<vector<TH2*>> *target=GetTarget(isData,charge);
    if(!target->size()||!target->at(0).size()){
      cout<<"[AddSetUpDown] no nominal hist"<<endl;
      exit(ENODATA);
    }
    TH2* central=target->at(0).at(0);
    vector<TH2*> hists;
    TString file=path(0,path.Index(':'));
    TString key=path(path.Index(':')+1,path.Length());
    if(!HasKey(file,key)){
      cout<<"[Efficiency::AddSetUpdown] no "<<key<<" in "<<file<<endl;
      exit(ENODATA);
    }
    TFile f(file);
    TH2* herror=(TH2*)f.Get(key);
    TH2* up=(TH2*)central->Clone();
    up->SetDirectory(NULL);
    for(int i=0;i<up->GetNcells();i++){
      up->SetBinContent(i,central->GetBinContent(i)+herror->GetBinContent(i));
    }
    hists.push_back(up);
    TH2* down=(TH2*)central->Clone();
    down->SetDirectory(NULL);
    for(int i=0;i<down->GetNcells();i++){
      down->SetBinContent(i,central->GetBinContent(i)-herror->GetBinContent(i));
    }
    hists.push_back(down);
    target->push_back(hists);
  }
}
void Efficiency::AddDataSetUpDown(TString path,int charge){
  AddSetUpDown(true,path,charge);
}
void Efficiency::AddSimSetUpDown(TString path,int charge){
  AddSetUpDown(false,path,charge);
}
const vector<vector<TH2*>>* Efficiency::GetTarget(bool isData,int charge) const{
  const vector<vector<TH2*>> *target=NULL;
  if(isData){
    if(charge>=0) target=&fDataPlus;
    else target=&fDataMinus;
  }else{
    if(charge>=0) target=&fSimPlus;
    else target=&fSimMinus;
  }
  return target;
}
vector<vector<TH2*>>* Efficiency::GetTarget(bool isData,int charge){
  return const_cast<vector<vector<TH2*>>*>(const_cast<const Efficiency*>(this)->GetTarget(isData,charge));
}
double Efficiency::GetEfficiency(bool isData,double eta,double pt,int charge,int set,int mem) const{
  const vector<vector<TH2*>> *target=GetTarget(isData,charge);
  int nset=target->size();
  if(nset<=set){
    cout<<"[Efficiency::GetEfficiency] no set "<<set<<endl;
    exit(ENODATA);
  }
  int nmem=target->at(set).size();
  if(nmem<=mem){
    cout<<"[Efficiency::GetEfficiency] no member "<<mem<<" for set "<<set<<endl;
    exit(ENODATA);
  }
  TH2* hist=target->at(set).at(mem);
  double xmin=hist->GetXaxis()->GetXmin();
  double xmax=hist->GetXaxis()->GetXmax();
  double ymin=hist->GetYaxis()->GetXmin();
  double ymax=hist->GetYaxis()->GetXmax();
  if(xmin>=0) eta=fabs(eta);
  if(eta<xmin) eta=xmin+0.001;
  if(eta>=xmax) eta=xmax-0.001;
  if(ymin>=0) pt=fabs(pt);
  if(pt<ymin) pt=ymin+0.001;
  if(pt>=ymax) pt=ymax-0.001;
  return hist->GetBinContent(hist->FindBin(eta,pt));
}
double Efficiency::GetDataEfficiency(double eta,double pt,int charge,int set,int mem) const{
  return GetEfficiency(true,eta,pt,charge,set,mem);
}
double Efficiency::GetSimEfficiency(double eta,double pt,int charge,int set,int mem) const{
  return GetEfficiency(false,eta,pt,charge,set,mem);
}
double Efficiency::GetEfficiencySF(double eta,double pt,int charge,int set,int mem) const{
  double data=GetEfficiency(true,eta,pt,charge,set,mem);
  double sim=GetEfficiency(false,eta,pt,charge,set,mem);
  if(sim==0) return 1.;
  return data/sim;
}
void Efficiency::Print(TString opt) const{
  cout<<"Data+: ";
  for(auto v:fDataPlus) cout<<v.size()<<" ";
  cout<<endl;
  if(opt.Contains("val")){
    for(int i=0;i<(int)fDataPlus.size();i++){
      for(int j=0;j<(int)fDataPlus[i].size();j++){
	cout<<GetDataEfficiency(2.4,200,+1,i,j)<<" ";
      }
      cout<<endl;
    }
  }
  cout<<"Data-: ";
  for(auto v:fDataMinus) cout<<v.size()<<" ";
  cout<<endl;
  if(opt.Contains("val")){
    for(int i=0;i<(int)fDataMinus.size();i++){
      for(int j=0;j<(int)fDataMinus[i].size();j++){
	cout<<GetDataEfficiency(2.4,200,-1,i,j)<<" ";
      }
      cout<<endl;
    }
  }
  cout<<"Sim+: ";
  for(auto v:fSimPlus) cout<<v.size()<<" ";
  cout<<endl;
  if(opt.Contains("val")){
    for(int i=0;i<(int)fSimPlus.size();i++){
      for(int j=0;j<(int)fSimPlus[i].size();j++){
	cout<<GetSimEfficiency(2.4,200,+1,i,j)<<" ";
      }
      cout<<endl;
    }
  }
  cout<<"Sim-: ";
  for(auto v:fSimMinus) cout<<v.size()<<" ";
  cout<<endl;
  if(opt.Contains("val")){
    for(int i=0;i<(int)fSimMinus.size();i++){
      for(int j=0;j<(int)fSimMinus[i].size();j++){
	cout<<GetSimEfficiency(2.4,200,-1,i,j)<<" ";
      }
      cout<<endl;
    }
  }
}
bool Efficiency::HasKey(TString path,TString key){
  if(!IsExists(path)){
    cout<<"[Efficiency::HasKey] no file "<<path<<endl;
    exit(ENODATA);
  }
  TFile f(path);
  TObject* obj=f.Get(key);
  if(obj){
    delete obj;
    return true;
  }
  return false;
}
bool Efficiency::IsExists(TString path){
  ifstream fcheck(path);
  return fcheck.good();
}
bool Efficiency::SameContents(TH1* hist1,TH1* hist2){
  if(!hist1) return false;
  if(!hist2) return false;
  if(hist1->GetNcells()!=hist2->GetNcells()) return false;
  for(int i=0;i<hist1->GetNcells();i++)
    if(hist1->GetBinContent(i)!=hist2->GetBinContent(i)) return false;
  return true;
}
EfficiencyTool::EfficiencyTool(TString path){
  if(path!="") Setup(path);
}
EfficiencyTool::~EfficiencyTool(){
  for(auto it:fEfficiencies){
    if(it.second) delete it.second;
  }
  fEfficiencies.clear();
}
void EfficiencyTool::Auto(TString key,TString path){
  Efficiency*& eff=fEfficiencies[key];
  if(!eff) eff=new Efficiency;
  int charge=0;
  if(IsPlus(path)) charge=+1;
  else if(IsMinus(path)) charge=-1;
  int nreplica=20;
  if(Efficiency::HasKey(path,"EGamma_EffData2D")){
    if(Efficiency::HasKey(path,"EGamma_EffData2D_stat")){
      eff->AddDataSetReplica(path+":EGamma_EffData2D",path+":EGamma_EffData2D_stat",nreplica,charge);
      eff->AddDataSet({path+":EGamma_EffData2D_altBkg"},charge);
      eff->AddDataSet({path+":EGamma_EffData2D_altSig"},charge);
      eff->AddDataSet({path+":EGamma_EffData2D"},charge);
      eff->AddDataSet({path+":EGamma_EffData2D"},charge);
      eff->AddSimSetReplica(path+":EGamma_EffMC2D",path+":EGamma_EffMC2D_stat",nreplica,charge);
      eff->AddSimSet({path+":EGamma_EffMC2D"},charge);
      eff->AddSimSet({path+":EGamma_EffMC2D"},charge);
      eff->AddSimSet({path+":EGamma_EffMC2D_altMC"},charge);
      eff->AddSimSet({path+":EGamma_EffMC2D_altTag"},charge);
    }else if(Efficiency::HasKey(path,"statData")){
      eff->AddDataSetReplica(path+":EGamma_EffData2D",path+":statData",nreplica,charge);
      eff->AddDataSetUpDown(path+":altBkgModel",charge);
      eff->AddDataSetUpDown(path+":altSignalModel",charge);
      eff->AddDataSet({path+":EGamma_EffData2D",path+":EGamma_EffData2D"},charge);
      eff->AddDataSet({path+":EGamma_EffData2D",path+":EGamma_EffData2D"},charge);
      eff->AddSimSetReplica(path+":EGamma_EffMC2D",path+":statMC",nreplica,charge);
      eff->AddSimSet({path+":EGamma_EffMC2D",path+":EGamma_EffMC2D"},charge);
      eff->AddSimSet({path+":EGamma_EffMC2D",path+":EGamma_EffMC2D"},charge);
      eff->AddSimSetUpDown(path+":altMCEff",charge);
      eff->AddSimSetUpDown(path+":altTagSelection",charge);
    }else{
      eff->AddDataSet({path+":EGamma_EffData2D"},charge);
      eff->AddSimSet({path+":EGamma_EffMC2D"},charge);
    } 
  }else if(Efficiency::HasKey(path,"muonEffi_data_eta_pt")){
    if(Efficiency::HasKey(path,"Systematics_data")){
      eff->AddDataSetReplica(path+":muonEffi_data_eta_pt",path+":muonEffi_data_eta_pt",nreplica,charge);
      eff->AddDataSet({path+":Systematics_data_massnarrow",path+":Systematics_data_massbroad"},charge);
      eff->AddDataSet({path+":Systematics_data_tagiso010",path+":Systematics_data_tagiso020"},charge);
      eff->AddSimSetReplica(path+":muonEffi_mc_eta_pt",path+":muonEffi_mc_eta_pt",nreplica,charge);
      eff->AddSimSet({path+":Systematics_mc_massnarrow",path+":Systematics_mc_massbroad"},charge);
      eff->AddSimSet({path+":Systematics_mc_tagiso010",path+":Systematics_mc_tagiso020"},charge);
      if(Efficiency::HasKey(path,"Systematics_data_altsig")){
	eff->AddDataSet({path+":Systematics_data_massbin50",path+":Systematics_data_massbin75"},charge);
	eff->AddDataSet({path+":Systematics_data_altsig"},charge);
	eff->AddSimSet({path+":Systematics_mc_massbin50",path+":Systematics_mc_massbin75"},charge);
	eff->AddSimSet({path+":Systematics_mc_altsig"},charge);
      }
    }else{
      eff->AddDataSet({path+":muonEffi_data_eta_pt"},charge);
      eff->AddSimSet({path+":muonEffi_mc_eta_pt"},charge);
    } 
  }
}
void EfficiencyTool::Auto(TString key,TString path1,TString path2){
  TString plus_path,minus_path;
  if(IsPlus(path1)&&IsMinus(path2)){
    Auto(key,path1);
    Auto(key,path2);
  }else if(IsPlus(path2)&&IsMinus(path1)){
    Auto(key,path2);
    Auto(key,path1);
  }else{
    cout<<"[EfficiencyTool::Auto] cannot determine charge "<<path1<<", "<<path2<<endl;
    exit(EXIT_FAILURE);
  }
}
void EfficiencyTool::Setup(TString path){
  if(!Efficiency::IsExists(path)){
    cout<<"[EfficiencyTool::Setup] no file "<<path<<endl;
    exit(ENODATA);
  }
  ifstream config(path);
  TString dirname=gSystem->DirName(path);
  string line;
  while(getline(config,line)){
    stringstream ss(line);
    vector<TString> words;
    {
      string word_temp;
      while(ss>>word_temp){
	TString word=word_temp;
	if(word.BeginsWith("#")) break;
	words.push_back(word);
      }
    }
    if(!words.size()) continue;
    if(words[0]=="auto"){
      if(words.size()==4){
	TString key=words[1];
	if(!words[2].BeginsWith("/")) words[2]=dirname+"/"+words[2];
	if(!words[3].BeginsWith("/")) words[3]=dirname+"/"+words[3];
	Auto(key,words[2],words[3]);
      }else if(words.size()==3){
	TString key=words[1];
	if(!words[2].BeginsWith("/")) words[2]=dirname+"/"+words[2];
	Auto(key,words[2]);
      }else{
	cout<<"[EfficiencyTool::Setup] wrong syntax words.size()="<<words.size()<<endl;
	cout<<line<<endl;
	exit(EXIT_FAILURE);
      }
    }else{
      cout<<"[EfficiencyTool::Setup] unknown command "<<words[0]<<endl;
      cout<<line<<endl;
      exit(EXIT_FAILURE);
    }
  }
  for(auto [key,eff]:fEfficiencies){
    cout<<"[EfficiencyTool::Setup] "<<key<<endl;
    eff->Print();
    //eff->Print("val");
  }
}
const Efficiency* EfficiencyTool::Get(TString key) const{
  const Efficiency* eff=NULL;
  if(fEfficiencies.find(key)!=fEfficiencies.end()) eff=fEfficiencies.find(key)->second;
  return eff;
}
double EfficiencyTool::GetDataEfficiency(TString key,double eta,double pt,int charge,int set,int mem) const{
  if(key==""||key=="Default") return 1.;
  const Efficiency* eff=Get(key);
  if(!eff){
    cout<<"[EfficiencyTool::GetDataEfficiency] no key "<<key<<endl;
    exit(ENODATA);
  }
  return eff->GetDataEfficiency(eta,pt,charge,set,mem);
}
double EfficiencyTool::GetDataEfficiency(TString key,const Lepton* lep,int set,int mem) const{
  if(key==""||key=="Default") return 1.;
  double eta=0;
  double pt=0;
  int charge=0;
  if(lep->InheritsFrom("Electron")){
    const Electron* el=(const Electron*)lep;
    eta=el->scEta();
    pt=el->UncorrPt();
    charge=el->Charge();
  }else if(lep->InheritsFrom("Muon")){
    const Muon* mu=(const Muon*)lep;
    eta=mu->Eta();
    pt=mu->MiniAODPt();    
    charge=mu->Charge();
  }
  return GetDataEfficiency(key,eta,pt,charge,set,mem);
}
double EfficiencyTool::GetSimEfficiency(TString key,double eta,double pt,int charge,int set,int mem) const{
  if(key==""||key=="Default") return 1.;
  const Efficiency* eff=Get(key);
  if(!eff){
    cout<<"[EfficiencyTool::GetSimEfficiency] no key "<<key<<endl;
    exit(ENODATA);
  }
  return eff->GetSimEfficiency(eta,pt,charge,set,mem);
}
double EfficiencyTool::GetSimEfficiency(TString key,const Lepton* lep,int set,int mem) const{
  if(key==""||key=="Default") return 1.;
  double eta=0;
  double pt=0;
  int charge=0;
  if(lep->InheritsFrom("Electron")){
    const Electron* el=(const Electron*)lep;
    eta=el->scEta();
    pt=el->UncorrPt();
    charge=el->Charge();
  }else if(lep->InheritsFrom("Muon")){
    const Muon* mu=(const Muon*)lep;
    eta=mu->Eta();
    pt=mu->MiniAODPt();    
    charge=mu->Charge();
  }
  return GetSimEfficiency(key,eta,pt,charge,set,mem);
}
double EfficiencyTool::GetEfficiencySF(TString key,double eta,double pt,int charge,int set,int mem) const{
  if(key==""||key=="Default") return 1.;
  const Efficiency* eff=Get(key);
  if(!eff){
    cout<<"[EfficiencyTool::GetEfficiencySF] no key "<<key<<endl;
    exit(ENODATA);
  }
  return eff->GetEfficiencySF(eta,pt,charge,set,mem);
}
double EfficiencyTool::GetEfficiencySF(TString key,const Lepton* lep,int set,int mem) const{
  if(key==""||key=="Default") return 1.;
  double eta=0;
  double pt=0;
  int charge=0;
  if(lep->InheritsFrom("Electron")){
    const Electron* el=(const Electron*)lep;
    eta=el->scEta();
    pt=el->UncorrPt();
    charge=el->Charge();
  }else if(lep->InheritsFrom("Muon")){
    const Muon* mu=(const Muon*)lep;
    eta=mu->Eta();
    pt=mu->MiniAODPt();    
    charge=mu->Charge();
  }
  return GetEfficiencySF(key,eta,pt,charge,set,mem);
}
vector<vector<double>> EfficiencyTool::GetStructure(TString key) const{
  vector<vector<double>> out;
  if(key==""||key=="Default") return out;
  const Efficiency* eff=Get(key);
  if(!eff){
    cout<<"[EfficiencyTool::GetStructure] no key "<<key<<endl;
    exit(ENODATA);
  }
  int nset=eff->fDataPlus.size();
  for(int i=0;i<nset;i++)
    out.push_back(vector<double>(eff->fDataPlus.at(i).size(),1.));
  return out;
}
bool EfficiencyTool::IsPlus(TString path){
  path.ToLower();
  if(path.Contains("plus")) return true;
  return false;
}
bool EfficiencyTool::IsMinus(TString path){
  path.ToLower();
  if(path.Contains("minus")) return true;
  return false;
}
