#ifndef ElectroWeakAnalysis_Aepcor_H
#define ElectroWeakAnalysis_Aepcor_H

#include <iostream>
#include <boost/math/special_functions/erf.hpp>

struct CrystalBallE{
    static const double pi;
    static const double sqrtPiOver2;
    static const double sqrt2;

    double m;
    double s;
    double aL;
    double aH;
    double nL;
    double nH;

    //precalculate some variables for faster evaluation
    double B[2];
    double E1;
    double D;
    double N;
    double NA[2];
    double Ns;
    double NC[2];
    double F[2];
    double G[2];
    double K[2];
    double cdfa[2];

    CrystalBallE():m(0),s(1),aL(1000),aH(1000),nL(1000),nH(1000){
	init();
    }

    void init(){
	for(int i:{0,1}) B[i] = getB(i);
	E1 = getE(1);
	D = getD();
	N = getN();
	for(int i:{0,1}) K[i] = getK(i);
	for(int i:{0,1}) NA[i] = N*getA(i);
	Ns = N*s;
        for(int i:{0,1}) NC[i] = Ns*getC(i);
        for(int i:{0,1}) F[i] = 1-fa(i)*fa(i)/n(i);
        for(int i:{0,1}) G[i] = s*n(i)/fa(i);
        cdfa[0] = cdf(m-aL*s);
        cdfa[1] = cdf(m+aH*s);
    }

     double n(int i, double k=1) const{return i==0? 1+(nL-1)/k : 1+(nH-1)/k;}
     double fa(int i) const{return i==0? std::abs(aL) : std::abs(aH);}
     double ex(int i) const{return exp(-fa(i)*fa(i)/2);}
     double getA(int i, double k=1) const{return pow(n(i,k)/fa(i),n(i,k))*ex(i);}
     double getC(int i, double k=1) const{return n(i,k)/fa(i)/(n(i,k)-1)*ex(i);}
     double getB(int i, double k=1) const{return n(i,k)/fa(i)-fa(i);}
     double getN(double k=1) const{return 1.0/s/(sqrtPiOver2*(erf(fa(0)/sqrt2) + erf(fa(1)/sqrt2)) + getC(0,k) + getC(1,k));}
     double getF(int i, double k=1) const {return 1-fa(i)*fa(i)/n(i,k);}
     double getG(int i, double k=1) const {return s*n(i,k)/fa(i);}
     double getD(double k=1) const{return sqrtPiOver2*erf(fa(0)/sqrt2)+getC(0,k);}
     double getK(int i, double k=1) const{return 1.0/(n(i,k)-1);}
     double getE(int i, double k=1) const{ return (sqrtPiOver2*(erf(fa(0)/sqrt2) + erf(fa(1)/sqrt2))+getC(0,k)+getC(1,k))/getC(i,k);}

    double pdf(double x, double ks=1) const{ 
	double d=(x-m)/(s*ks);
	if(d<-aL) return NA[0] / ks * pow(B[0]-d, -nL);
	if(d> aH) return NA[1] / ks * pow(B[1]+d, -nH);
	return           N / ks * exp(-d*d/2);
    }
    double cdf(double x, double ks=1) const{
	double d=(x-m)/(s*ks);
	if(d<-aL) return NC[0] / pow(F[0]-s*d/G[0], nL-1);
	if(d> aH) return NC[1] * (E1 - pow(F[1]+s*d/G[1], 1-nH) );
	return Ns * (D - sqrtPiOver2 * erf(-d/sqrt2));
    }
    double invcdf(double u, double ks=1) const{
	if(u<cdfa[0]) return m + ks*G[0]*(F[0] - pow(NC[0]/u, K[0]));
	if(u>cdfa[1]) return m - ks*G[1]*(F[1] - pow(E1-u/NC[1], -K[1]));
	return m - sqrt2 * s * ks* boost::math::erf_inv((D - u/Ns)/sqrtPiOver2);
    }

    double pdf(double x, double ks, double kt) const{ 
	double d=(x-m)/(s*ks);
	double tN = getN(kt);
        if(d<-aL) return tN * getA(0,kt) / ks * pow(getB(0,kt)-d, -n(0,kt));
        if(d> aH) return tN * getA(1,kt) / ks * pow(getB(1,kt)+d, -n(1,kt));
        return           tN / ks * exp(-d*d/2);
    }
    double cdf(double x, double ks, double kt) const{
	double d=(x-m)/(s*ks);
        if(d<-aL) return getN(kt)*s*getC(0,kt) / pow(getF(0,kt)-s*d/getG(0,kt), n(0,kt)-1);
        if(d> aH) return getN(kt)*s*getC(1,kt) * (getE(1,kt) - pow(getF(1,kt)+s*d/getG(1,kt), 1-n(1,kt)) );
        return getN(kt)* s * (getD(kt) - sqrtPiOver2 * erf(-d/sqrt2));
    }
    double invcdf(double u, double ks, double kt) const{
        if(u<cdf(m-s*aL)) return m + ks*getG(0,kt)*(getF(0,kt) - pow(getN(kt)*s*getC(0,kt)/u, getK(0,kt) ));
        if(u>cdf(m+s*aH)) return m - ks*getG(1,kt)*(getF(1,kt) - pow(getE(1,kt)-u/getN(kt)/s/getC(1,kt), -getK(1,kt)));
        return m - sqrt2 * s * ks* boost::math::erf_inv((getD(kt) - u/s/getN(kt))/sqrtPiOver2);
    }
};


struct Aepres{
    enum TYPE {MC, Data, Extra};

    struct ResParams{
	double eta; 
	std::vector<double> kRes; 
	std::vector<double> kTail; 
	std::vector<double> nTrk[2]; 
	std::vector<double> rsPar[5]; 
	std::vector<CrystalBallE> cb;
	ResParams():eta(0){}
    };

    int NETA;
    int NEXT;
    int NMIN;

    std::vector<ResParams> resol;

    Aepres();

    int etaBin(double x) const;
    int trkBin(double x, int h, TYPE T=MC) const;
    void reset();

    double rndm(int H, int F, double v) const;
    double Sigma(double pt, int H, int F) const;
    double kSpread(double gpt, double rpt, int ieta, int ir9, int iRun) const;
    double kSpread(double gpt, double rpt, int ieta) const;
    double nSmear(double pt, double eta, TYPE type, int n, double u) const;
    double uSmear(double pt, int ieta, TYPE type, int n, double u) const;
    double kSmear(double pt, double eta, TYPE type, double v, double u) const;
    double kSmear(double pt, double eta, TYPE type, double v, double u, int n) const;
    double kExtra(double pt, double eta, int nlayers, double u, double w) const;
    double kExtra(double pt, int ieta, int ir9, double u) const;
};

class Aepcor{

    private:
	enum TVAR{Default, Replica, Symhes};

	static const double MPHI; 

	int NETA;
	int NPHI; 
	int NRUN;

	int RETA;
	int NEXT;
	int NRUNRES;

	double DPHI;
	std::vector<double> etabin;
	std::vector<double> extbin;
	std::vector<int> runbin;
	std::vector<double> frunres;

	struct ShapeParams{
	    double M; 
	    double A;
	    double D;
	    ShapeParams(){reset();}
	    void reset(){M=1; A=0; D=0;}
	    double kCor(double pt) const{return M + A/pt - D*sqrt(pt);}
	};

	struct Aepone{
	    Aepres RR;
	    std::vector<std::vector<ShapeParams>> cShape[2];
	    std::vector<std::vector<double>> cPhi[2];
	    std::vector<std::vector<double>> cRun;
	};

	int nset;
	std::vector<int> nmem;
	std::vector<int> tvar;
	std::vector<std::vector<Aepone>> RC;
	template<typename T> int getBin(T x, const std::vector<T>& b) const;
	int phiBin(double phi) const;

    public:
	enum TYPE{MC, DT};
	Aepcor(); 
	Aepcor(std::string filename); 
	void init(std::string filename);
	void reset();

	const Aepres& getRes(int s=0, int m=0) const {return RC[s][m].RR;}
	bool empty() const{return RC.empty();}
	double getM(int T, int H, int F, int s=0, int m=0) const{return empty() ? 1.0 : RC[s][m].cShape[T][H][F].M;}
	double getA(int T, int H, int F, int s=0, int m=0) const{return empty() ? 0.0 : RC[s][m].cShape[T][H][F].A;}
	double getD(int T, int H, int F, int s=0, int m=0) const{return empty() ? 0.0 : RC[s][m].cShape[T][H][F].D;}
	double getP(int T, int H, int F, int s=0, int m=0) const{return empty() ? 1.0 : RC[s][m].cPhi[T][H][F];}
	double getR(int H, int F, int s=0, int m=0) const{return empty() ? 1.0 : RC[s][m].cRun[H][F];}

	double getK(int T, int H, int s=0, int m=0) const{return empty() ? 1.0 : RC[s][m].RR.resol[H].kRes[T];}
	double kGenSmear(double pt, double eta, double v, double u, Aepres::TYPE TT=Aepres::Data, int s=0, int m=0) const;
	double kScaleMC(double pt, double eta, double phi, double r9, int s=0, int m=0) const;

	double kScaleDTrun(double pt, double eta, double phi, double r9, int irun, int s=0, int m=0) const;
	double kScaleDT(double pt, double eta, double phi, double r9, int run, int s=0, int m=0) const;
	double kScaleDT4MC(double pt, double eta, double phi, double r9, double uRun, int s=0, int m=0) const;
	double kSpreadMC(double pt, double eta, double phi, double r9, double uRun, double gt, int s=0, int m=0) const;
	double kSmearMC(double pt, double eta, double phi, double r9, double u, int s=0, int m=0) const;
};

#endif
