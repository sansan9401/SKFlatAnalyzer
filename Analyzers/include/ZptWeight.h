#ifndef ZptWeight_h
#define ZptWeight_h

#include "SMPAnalyzerCore.h"

class ZptWeight : public SMPAnalyzerCore {

public:

  void executeEventWithChannelName(TString channelname);
  void executeEvent();
  void FillHists(Parameter& weights);

  ZptWeight();
  ~ZptWeight();

  static const int ptbinnum=52;
  static const int rapbinnum=6;
  double ptbin[ptbinnum+1]={0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,24,26,28,30,32,34,36,38,40,44,48,52,56,62,70,80,100,120,150,180,210,250,300,350,400,650};
  double rapbin[rapbinnum+1]={0.0,0.4,0.8,1.2,1.6,2.0,2.4};
  static const int genptbinnum=88;
  static const int genrapbinnum=24;
  double genptbin[genptbinnum+1]={0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,3.25,3.5,3.75,4,4.25,4.5,4.75,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18,18.5,19,19.5,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,42,44,46,48,50,52,54,56,59,62,66,70,75,80,90,100,200,400};
  double genrapbin[genrapbinnum+1]={0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4};
};



#endif

