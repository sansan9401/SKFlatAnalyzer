#ifndef Electron_h
#define Electron_h

#include "Lepton.h"
#include "TString.h"

class Electron: public Lepton{

public:

  Electron();
  ~Electron();

  void SetEnShift(double en_up, double en_down);
  inline double EnShift(int s) const {
    if(s==0) return 1.;
    else if(s>0) return j_En_up;
    else return j_En_down;
  }
  void SetResShift(double res_up, double res_down);
  inline double ResShift(int s) const {
    if(s==0) return 1.;
    else if(s>0) return j_Res_up;
    else return j_Res_down;
  }

  void SetSC(double sceta, double scphi, double sce);
  inline double scEta() const { return j_scEta; }
  inline double scPhi() const { return j_scPhi; }
  inline double scE() const { return j_scE; }
  void SetMVA(double mvaiso, double mvanoiso);
  inline double MVAIso() const { return j_mvaiso; }
  inline double MVANoIso() const { return j_mvanoiso; }

  void SetUncorrE(double une);
  inline double UncorrE() const { return j_EnergyUnCorr; }
  inline double UncorrPt() const { return Pt() * j_EnergyUnCorr/E(); }

  void SetPassConversionVeto(bool b);
  inline int PassConversionVeto() const { return j_passConversionVeto; }
  void SetNMissingHits(int n);
  inline int NMissingHits() const { return j_NMissingHits; };

  enum EtaRegion{
    IB, OB, GAP, EC
  };
  inline EtaRegion etaRegion() const {
    double sceta = fabs(scEta());
    if( sceta < 0.8 ) return IB;
    else if( sceta < 1.444 ) return OB;
    else if( sceta < 1.566 ) return GAP;
    else return EC;
  }
  
  void SetCutBasedIDVariables(
    double Full5x5_sigmaIetaIeta,
    double dEtaSeed,
    double dPhiIn,
    double HoverE,
    double InvEminusInvP,
    double e2x5OverE5x5,
    double e1x5OverE5x5,
    double trackIso,
    double dr03EcalRecHitSumEt,
    double dr03HcalDepth1TowerSumEt,
    double dr03HcalTowerSumEt,
    double dr03TkSumPt,
    double ecalPFClusterIso,
    double hcalPFClusterIso,
    int ecalDriven
  );
  inline double Full5x5_sigmaIetaIeta() const { return j_Full5x5_sigmaIetaIeta; }
  inline double dEtaSeed() const { return j_dEtaSeed; }
  inline double dPhiIn() const { return j_dPhiIn; }
  inline double HoverE() const { return j_HoverE; }
  inline double InvEminusInvP() const { return j_InvEminusInvP; }
  inline double e2x5OverE5x5() const { return j_e2x5OverE5x5; }
  inline double e1x5OverE5x5() const { return j_e1x5OverE5x5; }
  inline double TrkIso() const {return j_trkiso; }
  inline double dr03EcalRecHitSumEt() const { return j_dr03EcalRecHitSumEt; }
  inline double dr03HcalDepth1TowerSumEt() const { return j_dr03HcalDepth1TowerSumEt; }
  inline double dr03HcalTowerSumEt() const { return j_dr03HcalTowerSumEt; }
  inline double dr03TkSumPt() const { return j_dr03TkSumPt; }
  inline double ecalPFClusterIso() const { return j_ecalPFClusterIso; }
  inline double hcalPFClusterIso() const { return j_hcalPFClusterIso; }
  inline bool isEcalDriven() const { return j_isEcalDriven; }

  void SetIDBit(unsigned int idbit);

  void SetIDCutBit(vector<int> idcutbit);
  inline vector<int> IDCutBit() const { return j_IDCutBit; }

  static const int N_SELECTOR = 12;
  enum Selector {
    POG_CB_VETO = 1<< 0,
    POG_CB_LOOSE = 1<< 1,
    POG_CB_MEDIUM = 1<< 2,
    POG_CB_TIGHT = 1<< 3,
    POG_MVA_ISO_WP80 = 1<< 4,
    POG_MVA_ISO_WP90 = 1<< 5,
    POG_MVA_ISO_WPHZZ = 1<< 6,
    POG_MVA_ISO_WPLOOSE = 1<< 7,
    POG_MVA_NOISO_WP80 = 1<< 8,
    POG_MVA_NOISO_WP90 = 1<< 9,
    POG_MVA_NOISO_WPLOOSE = 1<< 10,
    POG_HEEP = 1<< 11
  };

  inline bool PassSelector( unsigned int s ) const { return (j_IDBit & s)==s; }

  inline bool passVetoID()   const {return PassSelector(POG_CB_VETO); }
  inline bool passLooseID()  const {return PassSelector(POG_CB_LOOSE); }
  inline bool passMediumID() const {return PassSelector(POG_CB_MEDIUM); }
  inline bool passTightID()  const {return PassSelector(POG_CB_TIGHT); }
  inline bool passMVAID_noIso_WP80() const {return PassSelector(POG_MVA_NOISO_WP80); }
  inline bool passMVAID_noIso_WP90() const {return PassSelector(POG_MVA_NOISO_WP90); }
  inline bool passMVAID_iso_WP80() const {return PassSelector(POG_MVA_ISO_WP80); }
  inline bool passMVAID_iso_WP90() const {return PassSelector(POG_MVA_ISO_WP90); }
  inline bool passHEEPID() const {return PassSelector(POG_HEEP); }
  bool passHEEP2018Prompt() const;

  bool Pass_SUSYMVAWP(TString wp) const;
  bool Pass_SUSYTight() const;
  bool Pass_SUSYLoose() const;

  void SetRelPFIso_Rho(double r);
  double EA();

  //==== ID
  bool PassID(TString ID) const;
  bool Pass_TESTID() const;

  bool Pass_CutBasedLooseNoIso() const;
  bool Pass_CutBasedVetoNoIso() const;
  bool Pass_CutBasedLoose() const;
  bool Pass_CutBasedVeto() const;
  void SetRho(double r);
  inline double Rho() const { return j_Rho; }
  void SetIsGsfCtfScPixChargeConsistent(bool b);
  inline bool IsGsfCtfScPixChargeConsistent() const { return j_isGsfCtfScPixChargeConsistent; }
  inline void SetR9(double r9) { j_r9=r9; }
  inline double R9() const { return j_r9; }
    
private:

  double j_En_up;
  double j_En_down;
  double j_Res_up;
  double j_Res_down;
  double j_scEta,j_scPhi, j_scE;
  double j_mvaiso, j_mvanoiso;
  double j_EnergyUnCorr;
  bool j_passConversionVeto;
  int j_NMissingHits;
  double j_Full5x5_sigmaIetaIeta, j_dEtaSeed, j_dPhiIn, j_HoverE, j_InvEminusInvP, j_e2x5OverE5x5, j_e1x5OverE5x5, j_trkiso, j_dr03EcalRecHitSumEt, j_dr03HcalDepth1TowerSumEt;
  double j_dr03HcalTowerSumEt, j_dr03TkSumPt, j_ecalPFClusterIso, j_hcalPFClusterIso;
  bool j_isEcalDriven;
  unsigned int j_IDBit;
  vector<int> j_IDCutBit;
  double j_RelPFIso_Rho;

  double j_Rho;
  int j_isGsfCtfScPixChargeConsistent;
  double j_r9;

  ClassDef(Electron,1)

};

#endif
