#ifndef EFFICIENCYTOOL_H
#define EFFICIENCYTOOL_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include "TSystem.h"
#include "TString.h"
#include "TH2.h"
#include "TFile.h"
#include "TRandom3.h"
#include "Electron.h"
#include "Muon.h"

using namespace std;

class Efficiency{
 public:
  vector<vector<TH2*>> fDataPlus;
  vector<vector<TH2*>> fDataMinus;
  vector<vector<TH2*>> fSimPlus;
  vector<vector<TH2*>> fSimMinus;
  
  Efficiency();
  ~Efficiency();

  void AddSet(bool isData,vector<TString> paths,int charge);
  void AddDataSet(vector<TString> paths,int charge);
  void AddSimSet(vector<TString> paths,int charge);
  void AddSetReplica(bool isData,TString nominal,TString stat,int nreplica,int charge);
  void AddDataSetReplica(TString nominal,TString stat,int nreplica,int charge);
  void AddSimSetReplica(TString nominal,TString stat,int nreplica,int charge);
  void AddSetUpDown(bool isData,TString path,int charge);
  void AddDataSetUpDown(TString path,int charge);
  void AddSimSetUpDown(TString path,int charge);
  
  const vector<vector<TH2*>>* GetTarget(bool isData,int charge) const;
  vector<vector<TH2*>>* GetTarget(bool isData,int charge);
  double GetEfficiency(bool isData,double eta,double pt,int charge,int set,int mem) const;
  double GetDataEfficiency(double eta,double pt,int charge,int set,int mem) const;
  double GetSimEfficiency(double eta,double pt,int charge,int set,int mem) const;
  double GetEfficiencySF(double eta,double pt,int charge,int set,int mem) const;

  void Print(TString opt="") const;
  static bool HasKey(TString path,TString key);
  static bool IsExists(TString path);
  static bool SameContents(TH1* hist1,TH1* hist2);
};
class EfficiencyTool{
 public:
  map<TString,Efficiency*> fEfficiencies;

  EfficiencyTool(TString path="");
  ~EfficiencyTool();
  
  void Auto(TString key,TString path);
  void Auto(TString key,TString path1,TString path2);
  void Setup(TString path);

  const Efficiency* Get(TString key) const;
  double GetDataEfficiency(TString key,double eta,double pt,int charge,int set,int mem) const;
  double GetDataEfficiency(TString key,const Lepton* lep,int set,int mem) const;
  double GetSimEfficiency(TString key,double eta,double pt,int charge,int set,int mem) const;
  double GetSimEfficiency(TString key,const Lepton* lep,int set,int mem) const;
  double GetEfficiencySF(TString key,double eta,double pt,int charge,int set,int mem) const;
  double GetEfficiencySF(TString key,const Lepton* lep,int set,int mem) const;
  vector<vector<double>> GetStructure(TString key) const;

  static bool IsPlus(TString path);
  static bool IsMinus(TString path);
};

#endif
