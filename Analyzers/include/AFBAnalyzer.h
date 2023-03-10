#ifndef AFBAnalyzer_h
#define AFBAnalyzer_h

#include "TKey.h"
#include "SMPAnalyzerCore.h"
#include "LHAPDF/Reweighting.h"

class AFBAnalyzer : public SMPAnalyzerCore {

public:
  virtual void test();
  virtual void initializeAnalyzer();
  virtual void executeEvent();
  virtual void executeEventGen();
  virtual void executeEventWithParameter(Parameter& p);
  virtual void executeEventWithParameter(Parameter&& p){Parameter pp=p;executeEventWithParameter(pp);}
  static int GetUnfoldBin(int nbin,const double* bins,double mass,double cost);
  virtual Parameter MakeParameter(TString key,TString option="");
  virtual bool PassSelection(Parameter& p);
  virtual void EvalWeights(Parameter& p);
  virtual void ResetRecoWeights(Parameter& p);
  virtual void FillHists(Parameter& p);

  AFBAnalyzer();
  ~AFBAnalyzer();

  //double GetCosThetaCS(const vector<Lepton*>& leps);
  virtual double GetCosThetaCS(const Particle *p0,const Particle *p1,int direction=0);
  virtual double GetCosThetaR(const Particle *l0,const Particle *l1,const Particle *j0,int direction=0);
  virtual double GetCosThetaT(const Particle *l0,const Particle *l1,const Particle *j0,int direction=0);
  virtual void FillHistsAFB(TString pre,TString hpre,TString suf,Particle* l0,Particle* l1,map<TString,double> map_weight);
  virtual void FillHardHists(TString pre,TString suf,const Gen& genparton0,const Gen& genparton1,const Gen& genhardl0,const Gen& genhardl1,const Gen& genhardj0,double w);
  //void FillGenAFBHists(TString pre,TString suf,const Gen& genl0,const Gen& genl1,const Gen& genphotons,double w);
  
  TString hardprefix;
  bool IsNominalRun=true;
  bool IsSkimmed=false;

  LHAPDF::PDF* PDFbase=NULL;
  LHAPDF::PDF* PDFnf4=NULL;
  
  static const int afb_mbinnum=40;
  static constexpr const double afb_mbin[afb_mbinnum+1]={52,56,60,65,70,74,77,80,82,84,86,88,89,90,91,92,93,94,96,98,100,103,106,110,115,120,130,140,150,175,200,240,280,340,400,500,600,700,800,1000,3000};
  static const int afb_ybinnum=12;
  static constexpr const double afb_ybin[afb_ybinnum+1]={-2.4,-2.0,-1.6,-1.2,-0.8,-0.4,0,0.4,0.8,1.2,1.6,2.0,2.4};
  static const int afb_ptbinnum=30;
  static constexpr const double afb_ptbin[afb_ptbinnum+1]={0,2,4,6,8,10,12,14,16,18,20,24,28,32,36,40,45,50,55,60,70,80,90,100,120,140,160,190,250,400,650};
  static const int afb_costbinnum=20;
  const double afb_costbin[afb_costbinnum+1]={-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};

  static const int grid_mbinnum=6;
  const double grid_mbin[grid_mbinnum+1]={52,60,80,100,150,400,3000};
  static const int grid_ybinnum=4;
  const double grid_ybin[grid_ybinnum+1]={-2.4,-1.2,0,1.2,2.4};
  static const int grid_ptbinnum=4;
  const double grid_ptbin[grid_ptbinnum+1]={0,20,50,100,650};

  static const int fine_mbinnum=38;
  const double fine_mbin[fine_mbinnum+1]={52,56,60,65,70,74,77,80,82,84,86,88,89,90,91,92,93,94,96,98,100,103,106,110,115,120,130,140,150,175,200,240,280,340,400,600,1000,2000,5000};
  static const int fine_ptbinnum=30;
  const double fine_ptbin[fine_ptbinnum+1]={0,2,4,6,8,10,12,14,16,18,20,24,28,32,36,40,45,50,55,60,70,80,90,100,120,140,160,190,250,500,1000};
  
  static const int lptbinnum=56;
  const double lptbin[lptbinnum+1]={0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,64,68,72,76,80,85,90,95,100,110,120,130,140,150,160,180,200,250,300,350,400,500,600,700,800,900,1000};
};



#endif

