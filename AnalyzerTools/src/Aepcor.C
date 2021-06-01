#ifndef ElectroWeakAnalysis_Aepcor
#define ElectroWeakAnalysis_Aepcor

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <exception>
#include "Aepcor.h"

const double CrystalBallE::pi = 3.14159;
const double CrystalBallE::sqrtPiOver2 = sqrt(CrystalBallE::pi/2.0);
const double CrystalBallE::sqrt2 = sqrt(2.0);

Aepres::Aepres(){
    reset();
}

void Aepres::reset(){
    NETA=0;
    NEXT=0;
    NMIN=0;
    std::vector<ResParams>().swap(resol);
}

int Aepres::etaBin(double eta) const{
    double abseta=fabs(eta);
    for(int i=0; i<NETA-1; ++i) if(abseta<resol[i+1].eta) return i;
    return NETA-1;
}

int Aepres::trkBin(double x, int h, TYPE T) const{
    for(int i=0; i<NEXT-1; ++i) if(x<resol[h].nTrk[T][i+1]) return i;
    return NEXT-1;
}

double Aepres::Sigma(double pt, int H, int F) const{
    const ResParams &rp = resol[H];
    double a = rp.rsPar[0][F];
    double b = rp.rsPar[1][F];
    double c = rp.rsPar[2][F];
    double e = sqrt(1.0 + b*b/pt + c*c/pt/pt);
    double d = rp.rsPar[3][F];
    double f = rp.rsPar[4][F];
    double t = 1.0 + d*pow(pt,f);
    return a*e*t/sqrt(e*e+t*t);
}

double Aepres::rndm(int H, int F, double w) const{
    const ResParams &rp = resol[H];
    return rp.nTrk[MC][F]+(rp.nTrk[MC][F+1]-rp.nTrk[MC][F])*w; 
}

double Aepres::kSpread(double gpt, double rpt, int ieta, int ir9, int iRun) const{
    double kold=rpt / gpt;
    double knew=1;
    const ResParams &rp = resol[ieta];
    double s = Sigma(gpt,ieta,ir9); //better in case of non-zero mean
    double u = 0;
    try{
	u = rp.cb[ir9].cdf((kold-1.0)/s, rp.kRes[MC], rp.kTail[MC]); 
	knew = 1.0 + s*rp.cb[ir9].invcdf(u, rp.kRes[1+iRun], rp.kTail[1+iRun]);
    }
    catch(std::exception e){
	std::cout << "Error: " << iRun << " " << gpt << " " << rpt << " " << s << " " << kold << " " << knew << " " << 
	    rp.kRes[MC] << " " << rp.kRes[1+iRun] << " " << rp.kTail[MC] << " " << rp.kTail[1+iRun] << " " << u << std::endl;
	return 1.0;
    }
    if(knew<0) return 1.0;
    return knew/kold;
}


double Aepres::kSpread(double gpt, double rpt, int ieta) const{
    const auto &k = resol[ieta].kRes;
    double x = rpt / gpt;
    return (1.0 + (x - 1.0)*k[Data]/k[MC]) / x;
}

double Aepres::nSmear(double pt, double eta, TYPE type, int n, double u) const{
    int H = etaBin(fabs(eta));
    int F = n>NMIN ? n-NMIN : 0;
    const ResParams &rp = resol[H];
    double x = rp.kRes[type] * Sigma(pt, H, F) * rp.cb[F].invcdf(u);
    return 1.0 + x;
}

double Aepres::uSmear(double pt, int H, TYPE type, int n, double u) const{
    int F = n>NMIN ? n-NMIN : 0;
    const ResParams &rp = resol[H];
    double x = rp.kRes[type] * Sigma(pt, H, F) * rp.cb[F].invcdf(u);
    return 1.0 + x;
}

double Aepres::kSmear(double pt, double eta, TYPE type, double v, double u) const{
    int H = etaBin(fabs(eta));
    int F = trkBin(v, H); 
    const ResParams &rp = resol[H];
    double x = rp.kRes[type] * Sigma(pt, H, F) * rp.cb[F].invcdf(u);
    return 1.0 + x;
}

double Aepres::kSmear(double pt, double eta, TYPE type, double w, double u, int n) const{
    int H = etaBin(fabs(eta));
    int F = n-NMIN;
    if(type==Data) F = trkBin(rndm(H, F, w), H, Data);
    const ResParams &rp = resol[H];
    double x = rp.kRes[type] * Sigma(pt, H, F) * rp.cb[F].invcdf(u);
    return 1.0 + x;
}

double Aepres::kExtra(double pt, double eta, int n, double u, double w) const{
    int H = etaBin(fabs(eta));
    int F = n>NMIN ? n-NMIN : 0;
    const ResParams &rp = resol[H];
    double v = rp.nTrk[MC][F]+(rp.nTrk[MC][F+1]-rp.nTrk[MC][F])*w;
    int D = trkBin(v, H, Data);
    double RD = rp.kRes[Data]*Sigma(pt, H, D);
    double RM = rp.kRes[MC]*Sigma(pt, H, F);
    double x = RD>RM ? sqrt(RD*RD-RM*RM)*rp.cb[F].invcdf(u) : 0;
    if(x<=-1) return 1.0;
    return 1.0 + x;
    return 1.0/(1.0 + x); 
}

double Aepres::kExtra(double pt, int ieta, int ir9, double u) const{
    const ResParams &rp = resol[ieta];
    double d = rp.kRes[Data];
    double m = rp.kRes[MC];
    double x = d>m ? sqrt(d*d-m*m) * Sigma(pt, ieta, ir9) * rp.cb[ir9].invcdf(u) : 0;
    if(x<=-1) return 1.0;
    return 1.0 + x;
}


Aepcor::Aepcor(){}

Aepcor::Aepcor(std::string filename){
    init(filename);
}

void Aepcor::reset(){
    NETA=0;
    NPHI=0;
    NRUN=0;
    std::vector<double>().swap(etabin);
    std::vector<double>().swap(extbin);
    std::vector<int>().swap(runbin);
    nset=0;
    std::vector<int>().swap(nmem);
    std::vector<std::vector<Aepone>>().swap(RC);
}


void Aepcor::init(std::string filename){
    std::ifstream in(filename.c_str());
    if(in.fail()) throw std::invalid_argument("Aepcor::init could not open file " + filename);

    int RMIN(0);
    RETA = 0;
    NEXT = 0;
    std::vector<double> BETA;

    std::string tag;
    int type, sys, mem, var, bin;	
    std::string s;
    while(std::getline(in, s)){
	std::stringstream ss(s); 
	std::string first4=s.substr(0,4);
	if(first4=="NSET"){
	    ss >> tag >> nset;
	    nmem.resize(nset);
	    tvar.resize(nset);
	    RC.resize(nset);
	}
	else if(first4=="NMEM") {
	    ss >> tag;
	    for(int i=0; i<nset; ++i) {
		ss >> nmem[i];
		RC[i].resize(nmem[i]);
	    }
	}
	else if(first4=="TVAR") {
	    ss >> tag;
	    for(int i=0; i<nset; ++i) ss >> tvar[i];
	}
	else if(first4=="RMIN") ss >> tag >> RMIN;
	else if(first4=="NEXT") {
	    ss >> tag >> NEXT;
	    extbin.resize(NEXT+1);
	    for(auto &h: extbin) ss >> h;
	}
	else if(first4=="RETA") {
	    ss >> tag >> RETA;
	    BETA.resize(RETA+1);
	    for(auto &h: BETA) ss >> h;

	}
	else if(first4=="CPHI") {
	    ss >> tag >> NPHI; 
	    DPHI=2*CrystalBallE::pi/NPHI;
	}
	else if(first4=="CETA")  {
	    ss >> tag >> NETA;
	    etabin.resize(NETA+1);
	    for(auto& h: etabin) ss >> h;
	}
	else if(first4=="KRUN")  {
	    ss >> tag >> NRUN;
	    runbin.resize(NRUN+1);
	    for(auto& h: runbin) ss >> h;
	}
	else if(first4=="FRUN")  {
	    ss >> tag >> NRUNRES;
	    frunres.resize(NRUNRES+1);
	    for(auto &h: frunres) ss >> h;
	}
	else{ 
	    ss >> sys >> mem >> tag;
	    auto &rc = RC[sys][mem]; 
	    rc.RR.NETA=RETA;
	    rc.RR.NEXT=NEXT;
	    rc.RR.NMIN=RMIN;
	    auto &resol = rc.RR.resol;
	    if(resol.empty()){
		resol.resize(RETA);
		for(size_t ir=0; ir<resol.size(); ++ir){
		    auto &r = resol[ir];
		    r.eta = BETA[ir];
		    r.kRes.resize(NRUNRES+1, 1.0);
		    r.kTail.resize(NRUNRES+1, 1.0);
		    r.cb.resize(NEXT);
		    for(auto i:{0,1})r.nTrk[i].resize(NEXT+1);
		    for(auto i:{0,1,2,3,4})r.rsPar[i].resize(NEXT);
		}
	    }

	    //init with 1's and 0's
	    auto &cs = rc.cShape;
	    for(TYPE T:{MC,DT}){
		if(cs[T].empty()){
		    cs[T].resize(RETA);
		    for(auto &i: cs[T]) {
			i.resize(NEXT);
			for(auto &k: i) k.reset();
		    }
		}
	    }
	    auto &cp = rc.cPhi;
	    for(TYPE T:{MC,DT}){
		if(cp[T].empty()){
		    cp[T].resize(NETA);
		    for(auto &i: cp[T]) i.resize(NPHI, 1.0);
		}
	    }

	    auto &ck = rc.cRun;
	    if(ck.empty()){
		ck.resize(RETA);
		for(auto& i: ck) i.resize(NRUN, 1.0);
	    }

	    if(tag=="R"){
		ss >> var >> bin; 
		for(int i=0; i<NEXT; ++i) {
		    switch(var){
			case 0: ss >> resol[bin].rsPar[var][i]; break;
			case 1: ss >> resol[bin].rsPar[var][i]; break;
			case 2: ss >> resol[bin].rsPar[var][i]; break; 
			case 3: ss >> resol[bin].rsPar[var][i]; break; 
			case 4: ss >> resol[bin].rsPar[var][i]; break; 
			case 5: ss >> resol[bin].cb[i].m; break; 
			case 6: ss >> resol[bin].cb[i].s; break; 
			case 7: ss >> resol[bin].cb[i].aL; break; 
			case 8: ss >> resol[bin].cb[i].nL; break; 
			case 9: ss >> resol[bin].cb[i].aH; break; 
			case 10: ss >> resol[bin].cb[i].nH; break; 
			default: break;
		    }
		}
	    }
	    else if(tag=="T") {
		ss >> type >> bin; 
		for(int i=0; i<NEXT+1; ++i) ss >> resol[bin].nTrk[type][i];
	    }
	    else if(tag=="F") {
		ss >> type; 
		for(int i=0; i<RETA; ++i) ss >> resol[i].kRes[type];

	    }
	    else if(tag=="L") {
		ss >> type; 
		for(int i=0; i<RETA; ++i) ss >> resol[i].kTail[type];
	    }
	    else if(tag=="S") {
		ss >> type >> var >> bin; 
		for(int i=0; i<NEXT; ++i){
		    auto &x = cs[type][bin][i];
		    if(var==0) { ss >> x.M; x.M = 1.0+x.M/100;}
		    else if(var==1){ ss >> x.A; }
		    else if(var==2){ ss >> x.D; x.D/=100; }
		}
	    }
	    else if(tag=="K") {
		ss >> bin; 
		for(int i=0; i<NRUN; ++i){
		    auto &x = ck[bin][i];
		    ss >> x;
		    x = 1.0 + x/100;
		}
	    }
	    else if(tag=="C") {
		ss >> type >>  bin; 
		for(int i=0; i<NPHI; ++i){
		    auto &x = cp[type][bin][i];
		    ss >> x;
		    x = 1.0 + x/100;
		}
	    }
	}
    }

    for(auto &rcs: RC)
       for(auto &rcm: rcs)
           for(auto &r: rcm.RR.resol)
               for(auto &i: r.cb) i.init();

    in.close();
}

const double Aepcor::MPHI=-CrystalBallE::pi;

template<typename T>
int Aepcor::getBin(T x, const std::vector<T>& b) const{
    if(x<=b[0]) return 0;
    auto it = upper_bound(b.begin(), b.end(), x);
    if(it==b.end()) return b.size()-2;
    return it - b.begin() -1;
}

int Aepcor::phiBin(double x) const{
    int ibin=(x-MPHI)/DPHI;
    if(ibin<0) return 0; 
    if(ibin>=NPHI) return NPHI-1;
    return ibin;
}

double Aepcor::kScaleDTrun(double pt, double eta, double phi, double r9, int irun, int s, int m) const{
    if(empty()) return 1.0;
    int ieta = getBin(eta, etabin);
    int iAbsEta = ieta>=NETA/2 ? ieta-NETA/2 : NETA/2-ieta-1; 
    int iphi = phiBin(phi);
    int ir9 = getBin(r9, extbin);
    const auto &r = RC[s][m]; 
    return r.cRun[iAbsEta][irun] * r.cShape[DT][iAbsEta][ir9].kCor(pt) * r.cPhi[DT][ieta][iphi];
}

double Aepcor::kScaleDT(double pt, double eta, double phi, double r9, int run, int s, int m) const{
    if(empty()) return 1.0;
    return kScaleDTrun(pt, eta, phi, r9, getBin(run, runbin), s, m); 
}

double Aepcor::kScaleDT4MC(double pt, double eta, double phi, double r9, double uRun, int s, int m) const{
    if(empty()) return 1.0;
    return kScaleDTrun(pt, eta, phi, r9, getBin(uRun, frunres), s, m); 
}


double Aepcor::kScaleMC(double pt, double eta, double phi, double r9, int s, int m) const{
    if(empty()) return 1.0;
    int ieta = getBin(eta, etabin);
    int iAbsEta = ieta>=NETA/2 ? ieta-NETA/2 : NETA/2-ieta-1; 
    int iphi = phiBin(phi);
    int ir9 = getBin(r9, extbin);
    const auto &r = RC[s][m]; 
    return r.cShape[DT][iAbsEta][ir9].kCor(pt) * r.cPhi[MC][ieta][iphi]; 
}

double Aepcor::kSpreadMC(double pt, double eta, double phi, double r9, double uRun, double gt, int s, int m) const{
    if(empty()) return 1.0;
    try{
	int ieta = getBin(eta, etabin);
	int iAbsEta = ieta>=NETA/2 ? ieta-NETA/2 : NETA/2-ieta-1; 
	int iphi = phiBin(phi);
	int ir9 = getBin(r9, extbin);
	const auto &r = RC[s][m];
	double k = r.cShape[MC][iAbsEta][ir9].kCor(pt) * r.cPhi[MC][ieta][iphi];
	int iRun = getBin(uRun, frunres);
	if(iRun>=NRUNRES) iRun=NRUNRES-1;
	return k * r.RR.kSpread(gt, k*pt, iAbsEta, ir9, iRun);
    }catch(std::exception e){
	std::cout << "input pars: " <<pt << " " << eta << " " << phi << " " << r9 << " " << uRun << " " << gt << " " << s << " " << m <<std::endl;
	throw e;
    }
}

double Aepcor::kSmearMC(double pt, double eta, double phi, double r9, double u, int s, int m) const{
    if(empty()) return 1.0;
    int ieta = getBin(eta, etabin);
    int iAbsEta = ieta>=NETA/2 ? ieta-NETA/2 : NETA/2-ieta-1; 
    int iphi = phiBin(phi);
    int ir9 = getBin(r9, extbin);
    const auto &r = RC[s][m];
    double k = r.cShape[MC][iAbsEta][ir9].kCor(pt) * r.cPhi[MC][ieta][iphi];
    return k * r.RR.kExtra(k*pt, iAbsEta, ir9, u);
}

double Aepcor::kGenSmear(double pt, double eta, double v, double u, Aepres::TYPE TT, int s, int m) const{
    if(empty()) return 1.0;
    return RC[s][m].RR.kSmear(pt, eta, TT, v, u);
}

#endif

