#include "MCCorrection.h"
#include "BTagEfficiencies/JetTagEfficiencies.C"

MCCorrection::MCCorrection() : 
IgnoreNoHist(false)
{

  histDir = TDirectoryHelper::GetTempDirectory("MCCorrection");

  genFinderDY = new GenFinderForDY();

}
vector<TString> MCCorrection::Split(TString s,TString del){
  TObjArray* array=s.Tokenize(del);
  vector<TString> out;
  for(const auto& obj:*array) out.push_back(((TObjString*)obj)->String());
  array->Delete();
  return out;
}
void MCCorrection::ReadHistograms(){

  TString datapath = getenv("DATA_DIR");

  TDirectory* origDir = gDirectory;

  //==== ID/Trigger
  TString IDpath = datapath+"/"+GetEra()+"/ID/";

  vector<TString> elhistmaps=Split(gSystem->GetFromPipe("find "+IDpath+"/Electron/ -name 'histmap*.txt' -type f"),"\n");
  for(const auto& elhistmap:elhistmaps){
    string elline;
    ifstream in(elhistmap);
    while(getline(in,elline)){
      std::istringstream is( elline );

      TString tstring_elline = elline;
      if(tstring_elline.Contains("#")) continue;

      TString a,b,c,d,e,f;
      is >> a; // ID,RECO
      is >> b; // Eff,SF
      is >> c; // <WPnames>
      is >> d; // <rootfilename>
      is >> e; // <histname>
      is >> f; // Class
      TFile *file = new TFile(IDpath+"/Electron/"+d);
      
      if(f=="TH2F"){
	histDir->cd();
	map_hist_Electron[a+"_"+b+"_"+c] = (TH2F *)file->Get(e)->Clone();
      }
      else if(f=="TGraphAsymmErrors"){
	histDir->cd();
	map_graph_Electron[a+"_"+b+"_"+c] = (TGraphAsymmErrors *)file->Get(e)->Clone();
      }
      else{
	cout << "[MCCorrection::MCCorrection] Wrong class type : " << elline << endl;
      }
      file->Close();
      delete file;
      origDir->cd();
    }
  }

  cout << "[MCCorrection::MCCorrection] map_hist_Electron :" << endl;
  for(std::map< TString, TH2F* >::iterator it=map_hist_Electron.begin(); it!=map_hist_Electron.end(); it++){
    cout << "[MCCorrection::MCCorrection] key = " << it->first << endl;
  }
  cout << "[MCCorrection::MCCorrection] map_graph_Electron :" << endl;
  for(std::map< TString, TGraphAsymmErrors* >::iterator it=map_graph_Electron.begin(); it!=map_graph_Electron.end(); it++){
    cout << "[MCCorrection::MCCorrection] key = " << it->first << endl;
  }


  vector<TString> muhistmaps=Split(gSystem->GetFromPipe("find "+IDpath+"/Muon/ -name 'histmap*.txt' -type f"),"\n");
  for(const auto& muhistmap:muhistmaps){
    string elline2;
    ifstream in2(muhistmap);
    while(getline(in2,elline2)){
      std::istringstream is( elline2 );

      TString tstring_elline = elline2;
      if(tstring_elline.Contains("#")) continue;
      
      TString a,b,c,d,e;
      is >> a; // ID,RERCO
      is >> b; // Eff,SF
      is >> c; // <WPnames>
      is >> d; // <rootfilename>
      is >> e; // <histname>
      TFile *file = new TFile(IDpath+"/Muon/"+d);
      histDir->cd();
      map_hist_Muon[a+"_"+b+"_"+c] = (TH2F *)file->Get(e)->Clone();
      file->Close();
      delete file;
      origDir->cd();
    }
  }

  cout << "[MCCorrection::MCCorrection] map_hist_Muon :" << endl;
  for(std::map< TString, TH2F* >::iterator it=map_hist_Muon.begin(); it!=map_hist_Muon.end(); it++){
    cout << "[MCCorrection::MCCorrection] key = " << it->first << endl;
  }


  // == Get Prefiring maps
  TString PrefirePath  = datapath+"/"+GetEra()+"/Prefire/";

  string elline3;
  ifstream in3(PrefirePath+"/histmap.txt");
  while(getline(in3,elline3)){
    std::istringstream is( elline3 );

    TString tstring_elline = elline3;
    if(tstring_elline.Contains("#")) continue;

    TString a,b,c;
    is >> a; // Jet, Photon
    is >> b; // <rootfilename>
    is >> c; // <histname>
    
    TFile *file = new TFile(PrefirePath+b);
    histDir->cd();
    map_hist_prefire[a + "_prefire"] = (TH2F *)file->Get(c)->Clone();
    file->Close();
    delete file;
    origDir->cd();
  }


  // == Get Pileup Reweight maps
  TString PUReweightPath = datapath+"/"+GetEra()+"/PileUp/";

  string elline4;
  ifstream  in4(PUReweightPath+"/histmap.txt");
  while(getline(in4,elline4)){
    std::istringstream is( elline4 );

    TString tstring_elline = elline4;
    if(tstring_elline.Contains("#")) continue;

    TString a,b,c;
    is >> a; // sample name
    is >> b; // syst
    is >> c; // rootfile name

    //if(DataYear == 2017 && a!=MCSample) continue;
    
    TFile *file = new TFile(PUReweightPath+c);
    if( (TH1D *)file->Get(a+"_"+b) ){
      histDir->cd();
      map_hist_pileup[a+"_"+b+"_pileup"] = (TH1D *)file->Get(a+"_"+b)->Clone();
    }
    else{
      cout << "[MCCorrection::ReadHistograms] No : " << a + "_" + b << endl;
    }
    file->Close();
    delete file;
    origDir->cd();
  }
/*
  cout << "[MCCorrection::MCCorrection] map_hist_pileup :" << endl;
  for(std::map< TString, TH1D* >::iterator it=map_hist_pileup.begin(); it!=map_hist_pileup.end(); it++){
    cout << it->first << endl;
  }
*/
  
  // == Get Official DY Pt reweight maps
  TString DYPtReweightPath = datapath+"/"+GetEra()+"/DYPtReweight/Zpt_weights_"+TString::Itoa(DataYear,10)+".root";
  TFile *file_DYPtReweightPath = new TFile(DYPtReweightPath);
  histDir->cd();
  hist_DYPtReweight_2D = (TH2D *)file_DYPtReweightPath->Get("zptmass_weights")->Clone();
  file_DYPtReweightPath->Close();
  delete file_DYPtReweightPath;
  origDir->cd();
}

MCCorrection::~MCCorrection(){

  delete genFinderDY;

}

void MCCorrection::SetMCSample(TString s){
  MCSample = s;
}
void MCCorrection::SetIsDATA(bool b){
  IsDATA = b;
}
void MCCorrection::SetEventInfo(int r, int l, int e){
  run = r;
  lumi = l;
  event = e;
}
void MCCorrection::SetIsFastSim(bool b){
  IsFastSim = b;
}

double MCCorrection::MuonReco_SF(TString key, double eta, double p, int sys){

  //cout << "[MCCorrection::MuonReco_SF] key = " << key << endl;
  //cout << "[MCCorrection::MuonReco_SF] eta = " << eta << ", p = " << p << endl;

  if(key=="Default") return 1.;

  double value = 1.;
  double error = 0.;

  eta = fabs(eta);

  if(key=="HighPtMuonRecoSF"){

    //==== XXX this histogram uses P not Pt    

    //==== boundaries
    if(p<50.) p = 50.;
    if(p>=3500.) p = 3499.;
    if(eta>=2.4) eta = 2.39;

  }

  //cout << "[MCCorrection::MuonReco_SF] corrected eta = " << eta << ", p = " << p << endl;
  //cout << "[MCCorrection::MuonReco_SF] histname = " << "RECO_SF_"+key << endl;

  TH2F *this_hist = map_hist_Muon["RECO_SF_"+key];
  if(!this_hist){
    if(IgnoreNoHist) return 1.;
    else{
      cerr << "[MCCorrection::MuonReco_SF] No "<<"RECO_SF_"+key<<endl;
      exit(ENODATA);
    }
  }

  int this_bin = this_hist->FindBin(p, eta);
  value = this_hist->GetBinContent(this_bin);
  error = this_hist->GetBinError(this_bin);

  //cout << "[MCCorrection::MuonReco_SF] --> value = " << value << "\t" << ", error = " << error << endl;

  return value+double(sys)*error;

}

double MCCorrection::MuonID_SF(TString ID, double eta, double pt, int sys){

  if(ID=="Default") return 1.;

  //cout << "[MCCorrection::MuonID_SF] ID = " << ID << endl;
  //cout << "[MCCorrection::MuonID_SF] eta = " << eta << ", pt = " << pt << endl;

  double value = 1.;
  double error = 0.;

  if(DataYear!=2016){
    eta = fabs(eta);
  }

  if(ID=="NUM_TightID_DEN_genTracks" || ID=="NUM_HighPtID_DEN_genTracks"){
    //==== boundaries
    if(pt<20.) pt = 20.;
    if(pt>=120.) pt = 119.;
    if(eta>=2.4) eta = 2.39;
    if(eta<-2.4) eta = -2.4;
  }

  TH2F *this_hist = map_hist_Muon["ID_SF_"+ID];
  if(!this_hist){
    if(IgnoreNoHist) return 1.;
    else{
      cerr << "[MCCorrection::MuonID_SF] No "<<"ID_SF_"+ID<<endl;
      exit(ENODATA);
    }
  }

  int this_bin(-999);

  if(DataYear==2016){
    this_bin = this_hist->FindBin(eta,pt);
  }
  else{
    this_bin = this_hist->FindBin(pt,eta);
  }

  value = this_hist->GetBinContent(this_bin);
  error = this_hist->GetBinError(this_bin);

  //cout << "[MCCorrection::MuonID_SF] value = " << value << endl;

  return value+double(sys)*error;

}

double MCCorrection::MuonISO_SF(TString ID, double eta, double pt, int sys){

  if(ID=="Default") return 1.;

  //cout << "[MCCorrection::MuonISO_SF] eta = " << eta << ", pt = " << pt << endl;

  double value = 1.;
  double error = 0.;

  if(DataYear!=2016){
    eta = fabs(eta);
  }

  if(ID=="NUM_TightRelIso_DEN_TightIDandIPCut" || ID=="NUM_LooseRelTkIso_DEN_HighPtIDandIPCut"){
    //==== boundaries
    if(pt<20.) pt = 20.;
    if(pt>=120.) pt = 119.;
    if(eta>=2.4) eta = 2.39;
    if(eta<-2.4) eta = -2.4;
  }

  TH2F *this_hist = map_hist_Muon["ISO_SF_"+ID];
  if(!this_hist){
    if(IgnoreNoHist) return 1.;
    else{
      cerr << "[MCCorrection::MuonISO_SF] No "<<"ISO_SF_"+ID<<endl;
      exit(ENODATA);
    }
  }

  int this_bin(-999);

  if(DataYear==2016){
    this_bin = this_hist->FindBin(eta,pt);
  }
  else{
    this_bin = this_hist->FindBin(pt,eta);
  }

  value = this_hist->GetBinContent(this_bin);
  error = this_hist->GetBinError(this_bin);

  //cout << "[MCCorrection::MuonISO_SF] value = " << value << endl;

  return value+double(sys)*error;

}

double MCCorrection::MuonTrigger_Eff(TString ID, TString trig, int DataOrMC, double eta, double pt, int sys){

  //cout << "[MCCorrection::MuonTrigger_Eff] Called" << endl;

  if(ID=="Default") return 1.;
  if(trig=="Default") return 1.;

  //cout << "[MCCorrection::MuonTrigger_Eff] DataYear = " << DataYear << endl;
  //cout << "[MCCorrection::MuonTrigger_Eff] ID = " << ID << "\t" << "trig = " << trig << endl;
  //cout << "[MCCorrection::MuonTrigger_Eff] DataOrMC = " << DataOrMC << endl;
  //cout << "[MCCorrection::MuonTrigger_Eff] eta = " << eta << ", pt = " << pt << endl;

  double value = 1.;
  double error = 0.;

  eta = fabs(eta);

  //==== 2016
  if(DataYear==2016){
    if(trig=="IsoMu24"){
      if(pt<26.) return 1.; //FIXME
      if(eta>=2.4) eta = 2.39;

      if(pt>500.) pt = 499.;
    }
    else if(trig=="Mu50"){
      if(pt<52.) return 1.; //FIXME
      if(eta>=2.4) eta = 2.39;

      if(pt>1000.) pt = 999.;
    }
    else{

    }
  }
  else if(DataYear==2017){
    if(trig=="IsoMu27"){
      //==== FIXME MiniAODPt Pt
      //==== FIXME 28.9918  29.0363
      //==== FIXME This event pass pt>29GeV cut, but MiniAOD pt < 29 GeV
      //==== FIXME So when I return 0., SF goes nan.. let's return 1 for now..
      if(pt<29.) return 1.; //FIXME
      if(eta>=2.4) eta = 2.39;

      if(pt>1200.) pt = 1199.;
    }
    else if(trig=="Mu50"){
      if(pt<52.) return 1.; //FIXME
      if(eta>=2.4) eta = 2.39;

      if(pt>1000.) pt = 999.;
    }
    else{

    }
  }
  else if(DataYear==2018){
    if(trig=="IsoMu24"){
      if(pt<26.) return 1.; //FIXME
      if(eta>=2.4) eta = 2.39;

      if(pt>1200.) pt = 1199.;
    }
    else if(trig=="Mu50"){
      if(pt<52.) return 1.; //FIXME
      if(eta>=2.4) eta = 2.39;

      if(pt>1000.) pt = 999.;
    }
    else{

    }
  }
  else{
    cerr << "[MCCorrection::MuonTrigger_Eff] Wrong year : " << DataYear << endl;
    exit(ENODATA);
  }

  TString histkey = "Trigger_Eff_DATA_"+trig+"_"+ID;
  if(DataOrMC==1) histkey = "Trigger_Eff_MC_"+trig+"_"+ID;
  //cout << "[MCCorrection::MuonTrigger_Eff] histkey = " << histkey << endl;
  TH2F *this_hist = map_hist_Muon[histkey];
  if(!this_hist){
    if(IgnoreNoHist) return 1.;
    else{
      cerr << "[MCCorrection::MuonTrigger_Eff] No "<<histkey<<endl;
      exit(ENODATA);
    }
  }

  int this_bin = this_hist->FindBin(pt,eta);

  value = this_hist->GetBinContent(this_bin);
  error = this_hist->GetBinError(this_bin);

  //cout << "[MCCorrection::MuonTrigger_Eff] value = " << value << endl;

  return value+double(sys)*error;


}

double MCCorrection::MuonTrigger_SF(TString ID, TString trig, const std::vector<Muon>& muons, int sys){

  if(ID=="Default") return 1.;
  if(trig=="Default") return 1.;

  double value = 1.;

  if(trig=="IsoMu24" || trig=="IsoMu27" || trig=="Mu50"){

    double eff_DATA = 1.;
    double eff_MC = 1.;

    for(unsigned int i=0; i<muons.size(); i++){
      eff_DATA *= ( 1.-MuonTrigger_Eff(ID, trig, 0, muons.at(i).Eta(), muons.at(i).MiniAODPt(), sys) );
      eff_MC   *= ( 1.-MuonTrigger_Eff(ID, trig, 1, muons.at(i).Eta(), muons.at(i).MiniAODPt(), -sys) );
    }

    eff_DATA = 1.-eff_DATA;
    eff_MC = 1.-eff_MC;

    value = eff_DATA/eff_MC;

/*
    if(eff_DATA==0||eff_MC==0){
      cout << "==== Zero Trigger Eff ====" << endl;
      for(unsigned int i=0;i<muons.size();i++){
        cout << muons.at(i).MiniAODPt() << "\t" << muons.at(i).Pt() << endl;
      }
    }
*/

  }

  return value;

}

double MCCorrection::MuonTrigger_SF(TString ID, TString trig, const std::vector<Muon *>& muons, int sys){

  std::vector<Muon> muvec;
  for(unsigned int i=0; i<muons.size(); i++){
    Muon this_muon = *(muons.at(i));
    muvec.push_back( this_muon );
  }

  return MuonTrigger_SF(ID, trig, muvec, sys);

}

double MCCorrection::ElectronID_SF(TString ID, double sceta, double pt, int sys){

  if(ID=="Default") return 1.;

  double value = 1.;
  double error = 0.;

  if(pt<10.) pt = 10.;
  if(pt>=500.) pt = 499.;
  if(sceta>=2.5) sceta = 2.49;
  if(sceta<-2.5) sceta = -2.5;

  if( ID.Contains("HEEP") ){

    //==== https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaRunIIRecommendations#HEEPV7_0
    //==== summary: a simple robust ID designed to be safe for high electrons.
    //==== The Et evolution of this ID must be well described in the MC therefore this ID is designed so its scale factor is flat vs Et.
    //==== As a result the HEEP differs in that it provides just a single number for the barrel and a single number for the endcap.
    //==== * note there almost certainly will have to be a retune for 2018 due to HCAL data/MC disagreements
    //==== * 2018 prompt: expected Dec 2018

    bool IsBarrel = fabs(sceta) < 1.479;
    double this_SF(1.);
    double this_SF_staterr(0.); // absolute value
    double this_SF_systerr(0.); // absolute value
    double this_SF_err(0.);

    if(DataYear==2016){
      this_SF         = (IsBarrel ? 0.971 : 0.983);
      this_SF_staterr = (IsBarrel ? 0.001 : 0.001);

      if(IsBarrel) this_SF_systerr = (pt<90. ? 0.01 : min(1.+(pt-90.)*0.0022,3.)*0.01) * this_SF;
      else         this_SF_systerr = (pt<90. ? 0.01 : min(1.+(pt-90.)*0.0143,4.)*0.01) * this_SF;

      this_SF_err = sqrt(this_SF_staterr*this_SF_staterr+this_SF_systerr*this_SF_systerr);
    }
    else if(DataYear==2017){ 
      this_SF         = (IsBarrel ? 0.967 : 0.973);
      this_SF_staterr = (IsBarrel ? 0.001 : 0.002);
      
      if(IsBarrel) this_SF_systerr = (pt<90. ? 0.01 : min(1.+(pt-90.)*0.0022,3.)*0.01) * this_SF;
      else         this_SF_systerr = (pt<90. ? 0.02 : min(1.+(pt-90.)*0.0143,5.)*0.01) * this_SF;
      
      this_SF_err = sqrt(this_SF_staterr*this_SF_staterr+this_SF_systerr*this_SF_systerr);
    }
    else if(DataYear==2018){

      //==== https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaRunIIRecommendations#HEEPv7_0_2018Prompt

      this_SF         = (IsBarrel ? 0.969 : 0.984);
      this_SF_staterr = (IsBarrel ? 0.000 : 0.001);

      if(IsBarrel) this_SF_systerr = (pt<90. ? 0.01 : min(1.+(pt-90.)*0.0022,3.)*0.01) * this_SF;
      else         this_SF_systerr = (pt<90. ? 0.02 : min(1.+(pt-90.)*0.0143,5.)*0.01) * this_SF;

      this_SF_err = sqrt(this_SF_staterr*this_SF_staterr+this_SF_systerr*this_SF_systerr);
    }
    else{
      cerr << "[MCCorrection::ElectronID_SF] (Hist) Wrong year "<< DataYear << endl;
      exit(ENODATA);
    }

    return this_SF+double(sys)*this_SF_err;

  }
  else{

    TH2F *this_hist = map_hist_Electron["ID_SF_"+ID];
    if(!this_hist){
      if(IgnoreNoHist) return 1.;
      else{
        cerr << "[MCCorrection::ElectronID_SF] (Hist) No "<<"ID_SF_"+ID<<endl;
        exit(ENODATA);
      }
    }

    int this_bin = this_hist->FindBin(sceta,pt);
    value = this_hist->GetBinContent(this_bin);
    error = this_hist->GetBinError(this_bin);

    return value+double(sys)*error;

  }

}

double MCCorrection::ElectronReco_SF(double sceta, double pt, int sys){

  double value = 1.;
  double error = 0.;

  TString ptrange = "ptgt20";
  if(pt<20.) ptrange = "ptlt20";

  if(pt<10.) pt = 10.;
  if(pt>=500.) pt = 499.;
  if(sceta>=2.5) sceta = 2.49;
  if(sceta<-2.5) sceta = -2.5;

  TH2F *this_hist = map_hist_Electron["RECO_SF_"+ptrange];
  if(!this_hist){
    if(IgnoreNoHist) return 1.;
    else{
      cerr << "[MCCorrection::ElectronReco_SF] No "<<"RECO_SF_"+ptrange<<endl;
      exit(ENODATA);
    }
  }

  //cout << "[MCCorrection::ElectronReco_SF] " << this_hist->GetBinContent(1,1) << endl;

  int this_bin = this_hist->FindBin(sceta,pt);
  value = this_hist->GetBinContent(this_bin);
  error = this_hist->GetBinError(this_bin);

  return value+double(sys)*error;

}

double MCCorrection::ElectronTrigger_Eff(TString ID, TString trig, int DataOrMC, double sceta, double pt, int sys){

  if(ID=="Default") return 1.;
  if(trig=="Default") return 1.;

  //cout << "[MCCorrection::ElectronTrigger_Eff] ID = " << ID << "\t" << "trig = " << trig << endl;
  //cout << "[MCCorrection::ElectronTrigger_Eff] DataOrMC = " << DataOrMC << endl;
  //cout << "[MCCorrection::ElectronTrigger_Eff] sceta = " << sceta << ", pt = " << pt << endl;

  double value = 1.;
  double error = 0.;

  if(trig=="WREGammaTrigger"){

    if(pt<50.) pt = 50.;
    if(pt>=500.) pt = 499.;

    if(sceta<-2.5) sceta = -2.5;
    if(sceta>=2.5) sceta = 2.49;

    TString etaregion = "Barrel";
    if(fabs(sceta) > 1.566) etaregion = "EndCap";

    TString histkey = "Trigger_Eff_DATA_"+trig+"_"+ID+"_"+etaregion;
    if(DataOrMC==1) histkey = "Trigger_Eff_MC_"+trig+"_"+ID+"_"+etaregion;
    //cout << "[MCCorrection::ElectronTrigger_Eff] histkey = " << histkey << endl;
    TH2F *this_hist = map_hist_Electron[histkey];
    if(!this_hist){
      if(IgnoreNoHist) return 1.;
      else{
        cerr << "[MCCorrection::ElectronTrigger_Eff] No "<<histkey<<endl;
        exit(ENODATA);
      }
    }

    int this_bin = this_hist->FindBin(sceta, pt);

    value = this_hist->GetBinContent(this_bin);
    error = this_hist->GetBinError(this_bin);

  }

  //cout << "[MCCorrection::ElectronTrigger_Eff] value = " << value << endl;

  return value+double(sys)*error;

}

double MCCorrection::ElectronTrigger_SF(TString ID, TString trig, const std::vector<Electron>& electrons, int sys){

  if(ID=="Default") return 1.;
  if(trig=="Default") return 1.;

  double value = 1.;

  if(trig=="WREGammaTrigger"){

    double eff_DATA = 1.;
    double eff_MC = 1.;

    for(unsigned int i=0; i<electrons.size(); i++){
      eff_DATA *= ( 1.-ElectronTrigger_Eff(ID, trig, 0, electrons.at(i).scEta(), electrons.at(i).Pt(), sys) );
      eff_MC   *= ( 1.-ElectronTrigger_Eff(ID, trig, 1, electrons.at(i).scEta(), electrons.at(i).Pt(), -sys) );
    }

    eff_DATA = 1.-eff_DATA;
    eff_MC = 1.-eff_MC;

    value = eff_DATA/eff_MC;
    if(IsFastSim) value = eff_DATA;


/*
    if(eff_DATA==0||eff_MC==0){
      cout << "==== Zero Trigger Eff ====" << endl;
      for(unsigned int i=0;i<electrons.size();i++){
        cout << electrons.at(i).Pt() << "\t" << electrons.at(i).Pt() << endl;
      }
    }
*/

  }

  return value;

}

double MCCorrection::ElectronTrigger_SF(TString ID, TString trig, const std::vector<Electron *>& electrons, int sys){

  std::vector<Electron> muvec;
  for(unsigned int i=0; i<electrons.size(); i++){
    Electron this_electron = *(electrons.at(i));
    muvec.push_back( this_electron );
  }

  return ElectronTrigger_SF(ID, trig, muvec, sys);

}

double MCCorrection::GetPrefireWeight(const std::vector<Photon>& photons, const std::vector<Jet>& jets, int sys){

  double photon_weight = 1.;
  double jet_weight = 1.;
  
  TH2F *photon_hist = map_hist_prefire["Photon_prefire"];
  TH2F *jet_hist = map_hist_prefire["Jet_prefire"];

  
  for(unsigned int i_pho = 0; i_pho < photons.size(); i_pho++){
    Photon current_photon = photons.at(i_pho);
    double eta = current_photon.scEta();
    double pt = current_photon.Pt();
    
    int this_bin = photon_hist->FindBin(eta, pt);

    double this_eff = photon_hist->GetBinContent(this_bin);
    double this_efferr = photon_hist->GetBinError(this_bin);

    double current_weight = 1. - (this_eff + (double)sys * this_efferr );;
    photon_weight = photon_weight * current_weight;
  }
  
  for(unsigned int i_jet = 0; i_jet < jets.size(); i_jet++){
    Jet current_jet = jets.at(i_jet);
    double eta = current_jet.Eta();
    double pt = current_jet.Pt();
    
    int this_bin = jet_hist->FindBin(eta, pt);

    double this_eff = jet_hist->GetBinContent(this_bin);
    double this_efferr = jet_hist->GetBinError(this_bin);

    double current_weight = 1.- (this_eff + (double)sys * this_efferr);
    jet_weight = jet_weight * current_weight;
  }

  return photon_weight * jet_weight;;

}


double MCCorrection::GetPileUpWeightBySampleName(int N_pileup, int syst){
  
  int this_bin = N_pileup+1;
  if(N_pileup >= 100) this_bin=100;

  TString this_histname = MCSample;
  if(syst == 0){
    this_histname += "_central_pileup";
  }
  else if(syst == -1){
    this_histname += "_sig_down_pileup";
  }
  else if(syst == 1){
    this_histname += "_sig_up_pileup";
  }
  else{
    cerr << "[MCCorrection::GetPileUpWeightBySampleName] syst should be 0, -1, or +1" << endl;
    exit(ENODATA);
  }

  TH1D *this_hist = map_hist_pileup[this_histname];
  if(!this_hist){
    cerr << "[MCCorrection::GetPileUpWeightBySampleName] No " << this_histname << endl;
    exit(ENODATA);
  }

  return this_hist->GetBinContent(this_bin);

}

double MCCorrection::GetPileUpWeight(int N_pileup, int syst){

  TString this_histname = "MC_" + GetEra();
  if(syst == 0){
    this_histname += "_central_pileup";
  }
  else if(syst == -1){
    this_histname += "_sig_down_pileup";
  }
  else if(syst == 1){
    this_histname += "_sig_up_pileup";
  }
  else{
    cerr << "[MCCorrection::GetPileUpWeight] syst should be 0, -1, or +1" << endl;
    exit(ENODATA);
  }

  TH1D *this_hist = map_hist_pileup[this_histname];
  if(!this_hist){
    cerr << "[MCCorrection::GetPileUpWeight] No " << this_histname << endl;
    exit(ENODATA);
  }
  
  int this_bin = N_pileup+1;
  if(this_bin>this_hist->GetNbinsX()) this_bin=this_hist->GetNbinsX();

  //FIXME 
  double pu_weight=this_hist->GetBinContent(this_bin);
  if(pu_weight>5) pu_weight=5;
  return pu_weight;

}

double MCCorrection::GetTopPtReweight(const std::vector<Gen>& gens){
  //==== ref: https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting2017
  //==== Only top quarks in SM ttbar events must be reweighted, 
  //==== not single tops or tops from BSM production mechanisms.
  if(!MCSample.Contains("TT") || !MCSample.Contains("powheg")){
    return 1.;
  }
  //==== initialize with large number
  double toppt1=10000, toppt2=10000;
  bool found_top = false, found_atop = false;

  for(vector<Gen>::const_iterator genit=gens.begin(); genit!=gens.end(); genit++){
    
    if(genit->Status() == 22){
      if(genit->PID() == 6){
        toppt1= genit->Pt();
        found_top = true;
      }
      else if(genit->PID() == -6){
        toppt2= genit->Pt();
        found_atop = true;
      }
    }
    //==== after we found top pair, break the loop
    if(found_top && found_atop) break;
  }
  double pt_reweight = 1.;
  //==== if top pair is not found, return 1.
  //==== the measurement covers only the range pt(top)<=800GeV, otherwise, return 1.
  if(toppt1<=800 && toppt2 <=800){
    pt_reweight*=exp(0.0615-0.0005*toppt1);
    pt_reweight*=exp(0.0615-0.0005*toppt2);
    pt_reweight = sqrt(pt_reweight);
  }
  return pt_reweight;
}

double MCCorrection::GetOfficialDYReweight(const vector<Gen>& gens, int sys){

  genFinderDY->Find(gens);
  Particle genZ = genFinderDY->GenZ;

  double mZ = genZ.M();
  double ptZ = genZ.Pt();

  int bin_mass = hist_DYPtReweight_2D->GetXaxis()->FindBin(mZ);
  int bin_pt   = hist_DYPtReweight_2D->GetYaxis()->FindBin(ptZ);

  double value = hist_DYPtReweight_2D->GetBinContent( bin_mass, bin_pt );
  double error = hist_DYPtReweight_2D->GetBinError( bin_mass, bin_pt );

  return value+double(sys)*error;

}

void MCCorrection::SetJetTaggingParameters(std::vector<JetTagging::Parameters> v){
  jetTaggingPars = v;
}

void MCCorrection::SetupJetTagging(){

  if(IsDATA) return;

  TString datapath = getenv("DATA_DIR");
  TString btagpath = datapath+"/"+TString::Itoa(DataYear,10)+"/BTag/";

  std::map< string, BTagCalibration > tmp_map_BTagCalibration; //==== key = tagger+"_"+method

  for(unsigned int i=0; i<jetTaggingPars.size(); i++){
    //==== (DeepCSV,Medium,incl,comb

    cout << "[MCCorrection::SetJetTaggingParameters] Contructing BTagCalibrationReader with ";jetTaggingPars.at(i).Print();
    string this_tagger = JetTagging::TaggerToString(jetTaggingPars.at(i).j_Tagger);

    string this_wp = JetTagging::WPToString(jetTaggingPars.at(i).j_WP);
    BTagEntry::OperatingPoint op = BTagEntry::OP_LOOSE;
    if(this_wp=="Loose"){
      op = BTagEntry::OP_LOOSE;
    }
    else if(this_wp=="Medium"){
      op = BTagEntry::OP_MEDIUM;
    }
    else if(this_wp=="Tight"){
      op = BTagEntry::OP_TIGHT;
    }
    else{
      cerr << "[MCCorrection::ReadJetTaggingCVSs()] Wrong WP : " << this_wp << endl;
      exit(ENODATA);
    }
    //==== When using iterativefit method, use BTagEntry::OP_RESHAPING
    if(jetTaggingPars.at(i).j_MeasurmentType_Light==JetTagging::iterativefit ||
       jetTaggingPars.at(i).j_MeasurmentType_Heavy==JetTagging::iterativefit ){
      op = BTagEntry::OP_RESHAPING;
    }

    string this_mt_L = JetTagging::MeasurmentTypeToString(jetTaggingPars.at(i).j_MeasurmentType_Light);
    string this_mt_H = JetTagging::MeasurmentTypeToString(jetTaggingPars.at(i).j_MeasurmentType_Heavy);

    ifstream in(btagpath+"/cvsmap.txt");
    string btagline; // dummy
    while(getline(in,btagline)){
      std::istringstream is( btagline );

      TString tstring_btagline = btagline;
      if(tstring_btagline.Contains("#")) continue;

      int tmp_Year, tmp_Run_Start, tmp_Run_End;
      string tmp_tagger, tmp_pd, tmp_filename;

      is >> tmp_Year; // YEAR
      is >> tmp_tagger; // TAGGER
      is >> tmp_pd; // Period_dep
      is >> tmp_Run_Start; // Run_start
      is >> tmp_Run_End; // Run_end
      is >> tmp_filename; // csv file
      if( tmp_tagger != this_tagger ) continue;
      if( tmp_pd != "All") continue; // don't no period dep now

      //==== Get BTagCalibration objects for both heavy-method and light-method first
      //==== BTagCalibration() is where we read and addEntry from csv file, so it takes time
      //==== So let's not contruct it again if it already exists
      string this_key_L = tmp_tagger+"_"+this_mt_L;
      std::map< string, BTagCalibration >::const_iterator tmp_it_L = tmp_map_BTagCalibration.find( tmp_tagger+"_"+this_key_L );
      if(tmp_it_L==tmp_map_BTagCalibration.end()){
        tmp_map_BTagCalibration[tmp_tagger+"_"+this_key_L] = BTagCalibration(tmp_tagger, this_mt_L, btagpath.Data()+tmp_filename);
      }
      string this_key_H = tmp_tagger+"_"+this_mt_H;
      std::map< string, BTagCalibration >::const_iterator tmp_it_H = tmp_map_BTagCalibration.find( tmp_tagger+"_"+this_key_H );
      if(tmp_it_H==tmp_map_BTagCalibration.end()){
        tmp_map_BTagCalibration[tmp_tagger+"_"+this_key_H] = BTagCalibration(tmp_tagger, this_mt_H, btagpath.Data()+tmp_filename);
      }

      //==== Now, contructing BTagCalibrationReader obect

      std::vector<std::string> systvec_L = {"up", "down"};
      std::vector<std::string> systvec_C = {"up", "down"};
      std::vector<std::string> systvec_B = {"up", "down"};
      if(this_mt_L=="iterativefit"){
        systvec_L = {"up_hf","down_hf","up_jes","down_jes","up_lfstats1","down_lfstats1","up_lfstats2","down_lfstats2"};
        systvec_C = {"up_cferr1","down_cferr1","up_cferr2","down_cferr2"};
        systvec_B = {"up_hfstats1","down_hfstats1","up_hfstats2","down_hfstats2","up_lf","down_lf","up_jes","down_jes"};
      }

      //==== Load L
      map_BTagCalibrationReader[tmp_tagger+"_"+this_wp+"_L_"+this_mt_L] = new BTagCalibrationReader(op, "central", systvec_L);
      map_BTagCalibrationReader[tmp_tagger+"_"+this_wp+"_L_"+this_mt_L]->load( tmp_map_BTagCalibration[tmp_tagger+"_"+this_key_L], BTagEntry::FLAV_UDSG, this_mt_L);
      //==== Load C
      map_BTagCalibrationReader[tmp_tagger+"_"+this_wp+"_C_"+this_mt_H] = new BTagCalibrationReader(op, "central", systvec_C);
      map_BTagCalibrationReader[tmp_tagger+"_"+this_wp+"_C_"+this_mt_H]->load( tmp_map_BTagCalibration[tmp_tagger+"_"+this_key_H], BTagEntry::FLAV_C, this_mt_H);
      //==== Load B
      map_BTagCalibrationReader[tmp_tagger+"_"+this_wp+"_B_"+this_mt_H] = new BTagCalibrationReader(op, "central", systvec_B);
      map_BTagCalibrationReader[tmp_tagger+"_"+this_wp+"_B_"+this_mt_H]->load( tmp_map_BTagCalibration[tmp_tagger+"_"+this_key_H], BTagEntry::FLAV_B, this_mt_H);

    }


  } // END loop jetTaggingPars

  cout << "[MCCorrection::SetJetTaggingParameters] Printing all BTagCalibrationReader :" << endl;
  for(std::map< std::string, BTagCalibrationReader* >::iterator it=map_BTagCalibrationReader.begin(); it!=map_BTagCalibrationReader.end(); it++){
    cout << "[MCCorrection::SetJetTaggingParameters] key = " << it->first << endl;
  }


}

double MCCorrection::GetJetTaggingSF(JetTagging::Parameters jtp, int JetFlavor, double JetPt, double JetEta, double Jetdiscr, string Syst){

  if(IsDATA) return 1.;

  string this_tagger = JetTagging::TaggerToString( jtp.j_Tagger );
  string this_wp = JetTagging::WPToString( jtp.j_WP );
  string this_mt_L = JetTagging::MeasurmentTypeToString(jtp.j_MeasurmentType_Light);
  string this_mt_H = JetTagging::MeasurmentTypeToString(jtp.j_MeasurmentType_Heavy);

  string key = JetTagging::TaggerToString( jtp.j_Tagger )+"_"+this_wp;
  BTagEntry::JetFlavor jf = BTagEntry::FLAV_B;
  if(abs(JetFlavor)==5){
    key += "_B_"+this_mt_H;
    jf = BTagEntry::FLAV_B;
  }
  else if(abs(JetFlavor)==4){
    key += "_C_"+this_mt_H;
    jf = BTagEntry::FLAV_C;
  }
  else{
    key += "_L_"+this_mt_L;
    jf = BTagEntry::FLAV_UDSG;
  }

  std::map< std::string, BTagCalibrationReader* >::const_iterator it = map_BTagCalibrationReader.find(key);
  if(it== map_BTagCalibrationReader.end()){
    cerr << "[MCCorrection::GetJetTaggingSF] b tag SF map not found for key = " << key << endl;
    return 1.;
    exit(ENODATA);
  }

  double this_SF = it->second->eval_auto_bounds(Syst, jf, fabs(JetEta), JetPt, Jetdiscr);
  //cout << "[MCCorrection::GetJetTaggingSF] key = " << it->first << endl;
  //cout << "[MCCorrection::GetJetTaggingSF] Jet tagging parameter : ";jtp.Print();
  //printf("[MCCorrection::GetJetTaggingSF] (JetFlavor, JetPt, JetEta, Jetdiscr, Syst) = (%d, %f, %f, %f, %s)\n",JetFlavor,JetPt,JetEta,Jetdiscr,Syst.c_str());
  //cout << "--> SF = " << this_SF << endl;
  return this_SF;

}

double MCCorrection::GetJetTaggingCutValue(JetTagging::Tagger tagger, JetTagging::WP wp){

  if(DataYear==2016){ //2016 values: pre-Legacy (To be fixed)
    if(tagger==JetTagging::DeepCSV){
      if(wp==JetTagging::Loose)  return 0.2217;
      if(wp==JetTagging::Medium) return 0.6321;
      if(wp==JetTagging::Tight)  return 0.8953;
    }
    if(tagger==JetTagging::DeepJet){
      if(wp==JetTagging::Loose)  return 0.0614;
      if(wp==JetTagging::Medium) return 0.3093;
      if(wp==JetTagging::Tight)  return 0.7221;
    }
  }
  if(DataYear==2017){
    //if(tagger==JetTagging::CSVv2){
    //  if(wp==JetTagging::Loose)  return 0.5803;
    //  if(wp==JetTagging::Medium) return 0.8838;
    //  if(wp==JetTagging::Tight)  return 0.9693;
    //}
    if(tagger==JetTagging::DeepCSV){
      if(wp==JetTagging::Loose)  return 0.1355; // 0.1522;
      if(wp==JetTagging::Medium) return 0.4506; // 0.4941;
      if(wp==JetTagging::Tight)  return 0.7738; // 0.8001;
    }
    if(tagger==JetTagging::DeepJet){
      if(wp==JetTagging::Loose)  return 0.0532; // 0.0521;
      if(wp==JetTagging::Medium) return 0.3040; // 0.3033;
      if(wp==JetTagging::Tight)  return 0.7476; // 0.7489;
    }
  }
  if(DataYear==2018){
    if(tagger==JetTagging::DeepCSV){
      if(wp==JetTagging::Loose)  return 0.1208; // 0.1241;
      if(wp==JetTagging::Medium) return 0.4168; // 0.4184;
      if(wp==JetTagging::Tight)  return 0.7665; // 0.7527;
    }
    if(tagger==JetTagging::DeepJet){
      if(wp==JetTagging::Loose)  return 0.0490; // 0.0494;
      if(wp==JetTagging::Medium) return 0.2783; // 0.2770;
      if(wp==JetTagging::Tight)  return 0.7100; // 0.7264;
    }
  }

  cout << "[MCCorrection::GetJetTaggingCutValue] Wrong " << endl;
  cout << "[MCCorrection::GetJetTaggingCutValue] DataYear = " << DataYear << endl;
  cout << "[MCCorrection::GetJetTaggingCutValue] tagger = " << tagger << endl;
  cout << "[MCCorrection::GetJetTaggingCutValue] wp = " << wp << endl;
  exit(ENODATA);

  return 1;

}

double MCCorrection::GetBTaggingReweight_1a(const vector<Jet>& jets, JetTagging::Parameters jtp, string Syst){

  if(IsDATA) return 1.;

  double Prob_MC(1.), Prob_DATA(1.);
  for(unsigned int i=0; i<jets.size(); i++){
    double this_MC_Eff = GetMCJetTagEff(jtp.j_Tagger, jtp.j_WP, jets.at(i).hadronFlavour(), jets.at(i).Pt(), jets.at(i).Eta());
    double this_SF = GetJetTaggingSF(jtp,
                                     jets.at(i).hadronFlavour(),
                                     jets.at(i).Pt(),
                                     jets.at(i).Eta(),
                                     jets.at(i).GetTaggerResult(jtp.j_Tagger),
                                     Syst );
    double this_DATA_Eff = this_MC_Eff*this_SF;

    bool isTagged = jets.at(i).GetTaggerResult(jtp.j_Tagger) > GetJetTaggingCutValue(jtp.j_Tagger, jtp.j_WP);
    if(isTagged){
      Prob_MC *= this_MC_Eff;
      Prob_DATA *= this_DATA_Eff;
    }
    else{
      Prob_MC *= 1.-this_MC_Eff;
      Prob_DATA *= 1.-this_DATA_Eff;
    }
  }

  return Prob_DATA/Prob_MC;

}

double MCCorrection::GetBTaggingReweight_1d(const vector<Jet>& jets, JetTagging::Parameters jtp, string Syst){

  if(IsDATA) return 1.;

  if(jtp.j_MeasurmentType_Light!=JetTagging::iterativefit || 
     jtp.j_MeasurmentType_Heavy!=JetTagging::iterativefit){
    cout << "[MCCorrection::GetBTaggingReweight_1d] This method only works for iterativefit method" << endl;
    cout << "[MCCorrection::GetBTaggingReweight_1d] jtp.j_MeasurmentType_Light = " << jtp.j_MeasurmentType_Light << endl;
    cout << "[MCCorrection::GetBTaggingReweight_1d] jtp.j_MeasurmentType_Heavy = " << jtp.j_MeasurmentType_Heavy << endl;
    exit(ENODATA);
    return 1.;
  }

  double rew(1.);

  for(unsigned int i=0; i<jets.size(); i++){

    int abs_hadFlavour = abs(jets.at(i).hadronFlavour());
    TString tmp_Syst(Syst);

/*
    systvec_L = {"up_hf","down_hf","up_jes","down_jes","up_lfstats1","down_lfstats1","up_lfstats2","down_lfstats2"};
    systvec_C = {"up_cferr1","down_cferr1","up_cferr2","down_cferr2"};
    systvec_B = {"up_hfstats1","down_hfstats1","up_hfstats2","down_hfstats2","up_lf","down_lf","up_jes","down_jes"};
*/

    bool GoodSyst = false;
    if(abs_hadFlavour==5){
      if(      tmp_Syst.Contains(TRegexp("hfstats[1-2]$")) ) GoodSyst = true;
      else if( tmp_Syst.Contains(TRegexp("lf$"))           ) GoodSyst = true;
      else if( tmp_Syst.Contains(TRegexp("jes$"))          ) GoodSyst = true;
    }
    else if(abs_hadFlavour==4){
      if(      tmp_Syst.Contains(TRegexp("cferr[1-2]$"))   ) GoodSyst = true;
    }
    else{
      if(      tmp_Syst.Contains(TRegexp("hf$")) )           GoodSyst = true;
      else if( tmp_Syst.Contains(TRegexp("jes$")) )          GoodSyst = true;
      else if( tmp_Syst.Contains(TRegexp("lfstats[1-2]$")) ) GoodSyst = true;
    }

    if(!GoodSyst) tmp_Syst = "central";

    double this_SF = GetJetTaggingSF(jtp,
                                     jets.at(i).hadronFlavour(),
                                     jets.at(i).Pt(),
                                     jets.at(i).Eta(),
                                     jets.at(i).GetTaggerResult(jtp.j_Tagger),
                                     string(tmp_Syst) );
    rew *= this_SF;
  }

  return rew;

}

bool MCCorrection::IsBTagged_2a(JetTagging::Parameters jtp, const Jet& jet, string Syst){

  double this_discr = jet.GetTaggerResult(jtp.j_Tagger);
  double cutValue = GetJetTaggingCutValue(jtp.j_Tagger, jtp.j_WP);

  bool isBTagged = this_discr > cutValue;

  if(IsDATA) return isBTagged;

  //==== Set seed
  unsigned int runNum_uint  = static_cast <unsigned int> (run);
  unsigned int lumiNum_uint = static_cast <unsigned int> (lumi);
  unsigned int evNum_uint   = static_cast <unsigned int> (event);
  unsigned int jet0eta = uint32_t(fabs(jet.Eta())/0.01);
  int m_nomVar=1;
  std::uint32_t seed = jet0eta + m_nomVar + (lumiNum_uint<<10) + (runNum_uint<<20) + evNum_uint;

  TRandom3 rand_(seed);

  bool newBTag = isBTagged;

  //=== Get SF
  double Btag_SF =  GetJetTaggingSF(jtp,
                                    jet.hadronFlavour(),
                                    jet.Pt(),
                                    jet.Eta(),
                                    jet.GetTaggerResult(jtp.j_Tagger),
                                    Syst );


  if(Btag_SF == 1) return newBTag; //no correction needed

  //=== throw random number to apply correction
  float coin = rand_.Uniform(1.);
  if(Btag_SF > 1){
    //=== use this if SF>1

    if( !isBTagged ) {

      double Btag_eff = GetMCJetTagEff(jtp.j_Tagger, jtp.j_WP, jet.hadronFlavour(), jet.Pt(), jet.Eta());
      //=== fraction of jets that need to be upgraded
      float mistagPercent = (1.0 - Btag_SF) / (1.0 - (1./Btag_eff) );

      //=== upgrade to tagged
      if( coin < mistagPercent ) {newBTag = true;}
    }

  }else{
    //=== use this if SF<1
 
    //=== downgrade tagged to untagged
    if( isBTagged && coin > Btag_SF ) {newBTag = false;}

  }

  return newBTag;

}

