#include "LHE.h"

ClassImp(LHE)

LHE::LHE() : Particle() {
  j_index = -1;
  j_ID = 0;
  j_status = 0;
}

LHE::~LHE(){
}

void LHE::SetIndexIDStatus(int i, int p, int t){
  j_index = i;
  j_ID = p;
  if(j_ID==1||j_ID==3||j_ID==5) SetCharge(-1./3);
  else if(j_ID==-1||j_ID==-3||j_ID==-5) SetCharge(+1./3);
  else if(j_ID==2||j_ID==4||j_ID==6) SetCharge(+2./3);
  else if(j_ID==-2||j_ID==-4||j_ID==-6) SetCharge(-2./3);
  else if(j_ID==11||j_ID==13||j_ID==15) SetCharge(-1);
  else if(j_ID==-11||j_ID==-13||j_ID==-15) SetCharge(+1);
  else SetCharge(0);
  j_status = t;
}

void LHE::Print(){

  cout << "(Index, ID, Status, Pt, Eta, Phi, M) = " << Index() << "\t" << ID() << "\t" << Status() << "\t" << Pt() << "\t" << Eta() << "\t" << Phi() << "\t" << M() << endl;

}

