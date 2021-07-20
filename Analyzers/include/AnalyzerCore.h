#ifndef AnalyzerCore_h
#define AnalyzerCore_h

#include "TLorentzVector.h"
#include "TString.h"
#include "TMath.h"
#include "TH3.h"
#include <sstream>      

#include "SKFlatNtuple.h"
#include "Event.h"
#include "Particle.h"
#include "Gen.h"
#include "LHE.h"
#include "Lepton.h"
#include "Muon.h"
#include "Electron.h"
#include "Photon.h"
#include "JetTaggingParameters.h"
#include "Jet.h"
#include "FatJet.h"

#include "AnalyzerParameter.h"
#include "MCCorrection.h"
#include "PuppiSoftdropMassCorr.h"
#include "FakeBackgroundEstimator.h"
#include "CFBackgroundEstimator.h"
#include "GeneralizedEndpoint.h"
#include "GEScaleSyst.h"
#include "PDFReweight.h"

#define M_Z 91.1876
#define M_W 80.379

class AnalyzerCore: public SKFlatNtuple {

public:

  AnalyzerCore();
  ~AnalyzerCore();

  virtual void initializeAnalyzer(){

  };

  virtual void executeEvent(){

  };

  //==================
  //==== Get objects
  //==================

  //==== When GetAllMuons, we apply Rochester correciton
  //==== Then, Pt orderding can be changed
  //==== So, in the analysis code, do e.g.,
  // vector<Muon> muons = GetMuons(~~~);
  // std::sort(muons.begin(), muons.end(), PtComparing);
  //==== ** Recommend you to do the same for other objects (Electron, Jet, FatJet, ...) **
  inline static bool PtComparing(const Particle& p1, const Particle& p2){ return (p1.Pt() > p2.Pt()); }
  inline static bool PtComparingPtr(Particle* p1, Particle* p2){ return (p1->Pt() > p2->Pt()); }

  Event GetEvent();


  std::vector<Electron> GetAllElectrons();
  std::vector<Electron> GetElectrons(TString id, double ptmin, double fetamax);

  std::vector<Muon> GetAllMuons();
  std::vector<Muon> GetMuons(TString id, double ptmin, double fetamax);

  std::vector<Photon> GetAllPhotons();
  std::vector<Photon> GetPhotons(TString id, double ptmin, double fetamax);

  //==== If TightIso is set, it calculate ptcone
  //==== If UseMini is true, Lepton::RelIso() returns MiniRelIso
  std::vector<Lepton *> MakeLeptonPointerVector(const std::vector<Muon>& muons, double TightIso=-999, bool UseMini=false);
  std::vector<Lepton *> MakeLeptonPointerVector(const std::vector<Electron>& electrons, double TightIso=-999, bool UseMini=false);

  std::vector<Jet> GetAllJets();
  std::vector<Jet> GetJets(TString id, double ptmin, double fetamax);

  std::vector<FatJet> GetAllFatJets();
  std::vector<FatJet> GetFatJets(TString id, double ptmin, double fetamax);

  std::vector<Gen> GetGens();
  std::vector<LHE> GetLHEs();

  //===================================================
  //==== Get objects METHOD 2
  //==== Get AllObject in the begining, and apply cut
  //==================================================+

  std::vector<Electron> SelectElectrons(const std::vector<Electron>& electrons, TString id, double ptmin, double fetamax);

  std::vector<Muon> UseTunePMuon(const std::vector<Muon>& muons);
  std::vector<Muon> SelectMuons(const std::vector<Muon>& muons, TString id, double ptmin, double fetamax);

  std::vector<Jet> SelectJets(const std::vector<Jet>& jets, TString id, double ptmin, double fetamax);

  std::vector<FatJet> SelectFatJets(const std::vector<FatJet>& jets, TString id, double ptmin, double fetamax);

  //==================
  //==== Systematics
  //==================

  std::vector<Electron> ScaleElectrons(const std::vector<Electron>& electrons, int sys);
  std::vector<Electron> SmearElectrons(const std::vector<Electron>& electrons, int sys);

  std::vector<Muon> ScaleMuons(const std::vector<Muon>& muons, int sys);

  std::vector<Jet> ScaleJets(const std::vector<Jet>& jets, int sys);
  std::vector<Jet> SmearJets(const std::vector<Jet>& jets, int sys);

  std::vector<FatJet> ScaleFatJets(const std::vector<FatJet>& jets, int sys);
  std::vector<FatJet> SmearFatJets(const std::vector<FatJet>& jets, int sys);
  std::vector<FatJet> ScaleSDMassFatJets(const std::vector<FatJet>& jets, int sys);
  std::vector<FatJet> SmearSDMassFatJets(const std::vector<FatJet>& jets, int sys);

  //====================
  //==== Event Filters
  //====================

  bool PassMETFilter();

  //============
  //==== Tools
  //============

  //===== Estimators

  MCCorrection *mcCorr=NULL;
  PuppiSoftdropMassCorr *puppiCorr=NULL;
  FakeBackgroundEstimator *fakeEst=NULL;
  CFBackgroundEstimator *cfEst=NULL;
  void initializeAnalyzerTools();

  //==== Prefire
  double GetPrefireWeight(int sys);

  //==== PU Reweight
  double GetPileUpWeight(int N_pileup, int syst);

  //==== Muon GeneralizedEngpoint momentum scaling
  GeneralizedEndpoint *muonGE=NULL;
  GEScaleSyst *muonGEScaleSyst=NULL;

  //==== Using new PDF set
  PDFReweight *pdfReweight=NULL;
  double GetPDFWeight(LHAPDF::PDF* pdf_);
  //==== NewCentral/ProdCentral
  double GetPDFReweight();
  //==== NewErrorSet/ProdCentral
  double GetPDFReweight(int member);

  //================
  //==== Functions
  //================

  bool IsOnZ(double m, double width);
  double MT(TLorentzVector a, TLorentzVector b);
  double MT2(TLorentzVector a, TLorentzVector b, Particle METv, double METgap);
  double projectedMET(TLorentzVector a, TLorentzVector b, Particle METv);
  bool HasFlag(TString flag);
  std::vector<Muon> MuonWithoutGap(const std::vector<Muon>& muons);
  std::vector<Muon> MuonPromptOnly(const std::vector<Muon>& muons, const std::vector<Gen>& gens);
  std::vector<Muon> MuonUsePtCone(const std::vector<Muon>& muons);
  Muon MuonUsePtCone(const Muon& muon);
  Particle UpdateMET(const Particle& METv, const std::vector<Muon>& muons);
  std::vector<Muon> MuonApplyPtCut(const std::vector<Muon>& muons, double ptcut);
  std::vector<Electron> ElectronPromptOnly(const std::vector<Electron>& electrons, const std::vector<Gen>& gens);
  std::vector<Electron> ElectronUsePtCone(const std::vector<Electron>& electrons);
  Electron ElectronUsePtCone(const Electron& electron);
  std::vector<Electron> ElectronApplyPtCut(const std::vector<Electron>& electrons, double ptcut);
  std::vector<Jet> JetsAwayFromFatJet(const std::vector<Jet>& jets, const std::vector<FatJet>& fatjets, double mindr=1.0);
  std::vector<Jet> JetsVetoLeptonInside(const std::vector<Jet>& jets, const std::vector<Electron>& els, const std::vector<Muon>& mus, double dR=0.4);
  std::vector<FatJet> FatJetsVetoLeptonInside(const std::vector<FatJet>& jets, const std::vector<Electron>& els, const std::vector<Muon>& mus, double dR=0.8);
  std::vector<Jet> JetsAwayFromPhoton(const std::vector<Jet>& jets, const std::vector<Photon>& photons, double mindr);
  Particle AddFatJetAndLepton(const FatJet& fatjet, const Lepton& lep);

  //==== GenMatching

  void PrintGen(const std::vector<Gen>& gens);
  static Gen GetGenMatchedLepton(const Lepton& lep, const std::vector<Gen>& gens);
  static Gen GetGenMatchedPhoton(const Lepton& lep, const std::vector<Gen>& gens);
  static vector<int> TrackGenSelfHistory(const Gen& me, const std::vector<Gen>& gens);
  bool IsFromHadron(const Gen& me, const std::vector<Gen>& gens);
  int GetLeptonType(const Lepton& lep, const std::vector<Gen>& gens);
  int GetLeptonType_Public(int TruthIdx, const std::vector<Gen>& TruthColl);
  int GetGenPhotonType(const Gen& genph, const std::vector<Gen>& gens);
  static bool IsFinalPhotonSt23_Public(const std::vector<Gen>& TruthColl);
  int  GetPrElType_InSameSCRange_Public(int TruthIdx, const std::vector<Gen>& TruthColl);
  bool IsSignalPID(int pid);

  //==== Plotting

  std::map< TString, TH1D* > maphist_TH1D;
  std::map< TString, TH2D* > maphist_TH2D;
  std::map< TString, TH3D* > maphist_TH3D;

  TH1D* GetHist1D(TString histname);
  TH2D* GetHist2D(TString histname);
  TH3D* GetHist3D(TString histname);

  void FillHist(TString histname, double value, double weight, int n_bin, double x_min, double x_max);
  void FillHist(TString histname, double value, double weight, int n_bin, double *xbins);
  void FillHist(TString histname,
                double value_x, double value_y,
                double weight,
                int n_binx, double x_min, double x_max,
                int n_biny, double y_min, double y_max);
  void FillHist(TString histname,
                double value_x, double value_y,
                double weight,
                int n_binx, double *xbins,
                int n_biny, double *ybins);
  void FillHist(TString histname,
		double value_x, double value_y, double value_z,
		double weight,
		int n_binx, double x_min, double x_max,
		int n_biny, double y_min, double y_max,
		int n_binz, double z_min, double z_max);
  void FillHist(TString histname,
		double value_x, double value_y, double value_z,
		double weight,
		int n_binx, const double *xbins,
		int n_biny, const double *ybins,
		int n_binz, const double *zbins);

  //==== JSFillHist : 1D
  std::map< TString, std::map<TString, TH1D*> > JSmaphist_TH1D;
  TH1D* JSGetHist1D(TString suffix, TString histname);
  void JSFillHist(TString suffix, TString histname, double value, double weight, int n_bin, double x_min, double x_max);
  //==== JSFillHist : 2D
  std::map< TString, std::map<TString, TH2D*> > JSmaphist_TH2D;
  TH2D* JSGetHist2D(TString suffix, TString histname);
  void JSFillHist(TString suffix, TString histname,
                  double value_x, double value_y,
                  double weight,
                  int n_binx, double x_min, double x_max,
                  int n_biny, double y_min, double y_max);
  void JSFillHist(TString suffix, TString histname,
                  double value_x, double value_y,
                  double weight,
                  int n_binx, double *xbins,
                  int n_biny, double *ybins);

  virtual void WriteHist();

  //==== Quick Plotters
  void FillLeptonPlots(std::vector<Lepton *> leps, TString this_region, double weight);
  void FillJetPlots(std::vector<Jet> jets, std::vector<FatJet> fatjets, TString this_region, double weight);

  //==== Output rootfile
  void SwitchToTempDir();
  TFile *outfile=NULL;
  void SetOutfilePath(TString outname);

};



#endif

