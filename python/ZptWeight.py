import os
import scipy.integrate as integrate
from scipy.optimize import minimize
import ROOT as rt
from math import exp,sqrt,log
from array import array
rt.gROOT.LoadMacro("./Plotter/AFBPlotter.cc")

def Short(era):
    return era.replace("2016preVFP","2016a").replace("2016postVFP","2016b")

def GetHist(args,ientry,mrange=None,yrange=None,zptweight=False,gen=False):
    histname="myptgpt_nozptweight"
    if gen: histname="gen_mypt_nozptweight"
    option="noproject"
    hist=None
    for channel in args.channel:
        for era in args.era:
            if channel[0]=="m":
                h=args.mPlotter.GetHist(ientry,channel+Short(era)+"/"+histname,option)
            elif channel[0]=="e":
                h=args.ePlotter.GetHist(ientry,channel+Short(era)+"/"+histname,option)
            if hist is None: 
                hist=h
            else: 
                hist.Add(h)
                h.Delete()

    if mrange is None:
        mrange=[rt.ZptWeight.massbin[0],rt.ZptWeight.massbin[rt.ZptWeight.massbinnum]]
    if yrange is None:
        yrange=[rt.ZptWeight.ybin[0],rt.ZptWeight.ybin[rt.ZptWeight.ybinnum]]
    ixmin=hist.GetXaxis().FindBin(mrange[0])
    ixmax=hist.GetXaxis().FindBin(mrange[1]-1e-6)
    iymin=hist.GetYaxis().FindBin(yrange[0])
    iymax=hist.GetYaxis().FindBin(yrange[1]-1e-6)

    if zptweight:
        if hist.InheritsFrom("TH4D"):
            for ix in range(ixmin,ixmax+1):
                for iy in range(iymin,iymax+1):
                    for iu in range(hist.GetUaxis().GetNbins()+2):
                        m=hist.GetXaxis().GetBinCenter(ix)
                        y=hist.GetYaxis().GetBinCenter(iy)
                        pt=hist.GetUaxis().GetBinCenter(iu)
                        sf=GetZptWeight(args,m,y,pt)
                        for iz in range(hist.GetZaxis().GetNbins()+2):
                            hist.SetBinContent(ix,iy,iz,iu,sf*hist.GetBinContent(ix,iy,iz,iu))
                            hist.SetBinError(ix,iy,iz,iu,sf*hist.GetBinError(ix,iy,iz,iu))
        elif hist.InheritsFrom("TH3"):
            for ix in range(ixmin,ixmax+1):
                for iy in range(iymin,iymax+1):
                    for iz in range(hist.GetZaxis().GetNbins()+2):
                        m=hist.GetXaxis().GetBinCenter(ix)
                        y=hist.GetYaxis().GetBinCenter(iy)
                        pt=hist.GetZaxis().GetBinCenter(iz)
                        sf=GetZptWeight(args,m,y,pt)
                        hist.SetBinContent(ix,iy,iz,sf*hist.GetBinContent(ix,iy,iz))
                        hist.SetBinError(ix,iy,iz,sf*hist.GetBinError(ix,iy,iz))
                            
    hist1d=hist.ProjectionZ("pt",ixmin,ixmax,iymin,iymax)
    hist.Delete()
    return hist1d
    
def GetDataPt(args): return GetHist(args,0,"pt")
def GetDYPt(args): return GetHist(args,1,"pt")
def GetDYGenPt(args): return GetHist(args,1,"genpt")

def FuncZptWeight(xx,par):
    #npar=6
    kk=[0,20,50,200]
    x=xx[0]
    aa=[]
    bb=[]
    cc=[]
    dd=[]
    for i in range(len(kk)-1):
        if i==0:
            aa+=[par[0]]
            bb+=[par[1]]
            cc+=[par[2]]
            dd+=[par[3]]
        else:
            aa+=[aa[i-1]+bb[i-1]*(kk[i]-kk[i-1])+cc[i-1]*(kk[i]-kk[i-1])**2+dd[i-1]*(kk[i]-kk[i-1])**3]
            bb+=[bb[i-1]+2*cc[i-1]*(kk[i]-kk[i-1])+3*dd[i-1]*(kk[i]-kk[i-1])**2]
            #cc+=[par[2*i+2]]
            cc+=[2*cc[i-1]+6*dd[i-1]*(kk[i]-kk[i-1])]
            #dd+=[par[2*i+3]]
            dd+=[par[i+3]]
        if x<kk[i+1]:
            return aa[i]+bb[i]*(x-kk[i])+cc[i]*(x-kk[i])**2+dd[i]*(x-kk[i])**3
    return (aa[-1]+bb[-1]*(kk[-1]-kk[-2])+cc[-1]*(kk[-1]-kk[-2])**2+dd[-1]*(kk[-1]-kk[-2])**3)+(bb[-1]+2*cc[-1]*(kk[-1]-kk[-2])+3*dd[-1]*(kk[-1]-kk[-2])**2)*(x-kk[-1])

def SetupZptWeight(args):
    npar=6
    mbins=[52,80,90,100,400]
    ybins=[0.0,0.4,1.0,1.6,2.4]
    #mbins=[80,90,100]
    #ybins=[0,1.2,2.4]
    args.hbin=rt.TH2D("bins","bins",len(mbins)-1,array("d",mbins),len(ybins)-1,array("d",ybins))
    args.hbin.SetDirectory(0)
    args.zptweight=[None]*args.hbin.GetNcells()
    for i in range(1,args.hbin.GetNbinsX()+1):
        for j in range(1,args.hbin.GetNbinsY()+1):
            this_bin=args.hbin.GetBin(i,j)
            args.zptweight[this_bin]=rt.TF1("zptweight{}".format(this_bin),FuncZptWeight,0,650,npar)
            for k in range(npar):
                args.zptweight[this_bin].SetParameter(k,0)
            args.zptweight[this_bin].SetParameter(0,1)

def CheckZptWeight(args):
    cs=[]
    gs=[]
    ms=[60,90,300]
    ys=[0,1,2,4]
    for im in range(len(ms)):
        cs+=[rt.TCanvas()]
        for iy in range(len(ys)):
            npoints=400
            xpoints=[i for i in range(npoints)]
            ypoints=[GetZptWeight(args,ms[im],ys[iy],x) for x in xpoints]
            gs+=[rt.TGraph(npoints,array("d",xpoints),array("d",ypoints))]
            if iy==0:
                gs[-1].Draw("al")
                gs[-1].SetTitle("zptweight at m={}".format(ms[im]))
                gs[-1].GetYaxis().SetRangeUser(0.5,1.5)
            else:
                gs[-1].SetLineColor(iy+1)
                gs[-1].Draw("same l")
    cs+=[rt.TCanvas()]
    for im in range(len(ms)):
        npoints=400
        xpoints=[i for i in range(npoints)]
        ypoints=[GetZptWeight(args,ms[im],0,x) for x in xpoints]
        gs+=[rt.TGraph(npoints,array("d",xpoints),array("d",ypoints))]
        if im==0:
            gs[-1].Draw("al")
            gs[-1].SetTitle("zptweight at y=0")
            gs[-1].GetYaxis().SetRangeUser(0.5,1.5)
        else:
            gs[-1].SetLineColor(im+1)
            gs[-1].Draw("same l")
    cs+=[rt.TCanvas()]
    hframe=cs[-1].DrawFrame(0,0.5,400,1.5)
    for i in range(len(args.zptweight)):
        if args.zptweight[i] is None: continue
        args.zptweight[i].Draw("same")
    f=rt.TFile("out.root","recreate")
    for c in cs:
        c.Write()
    

def GetZptWeight(args,mass,rapidity,pt):
    x=mass
    y=rapidity
    xmin=args.hbin.GetXaxis().GetBinCenter(1)
    xmax=args.hbin.GetXaxis().GetBinCenter(args.hbin.GetXaxis().GetLast())
    if x<xmin:
        binx1=1
        binx2=2
    elif x>=xmax:
        binx1=args.hbin.GetXaxis().GetLast()-1
        binx2=args.hbin.GetXaxis().GetLast()
    else:
        binx=args.hbin.GetXaxis().FindBin(x)
        if x>=args.hbin.GetXaxis().GetBinCenter(binx):
            binx1=binx
            binx2=binx+1
        else:
            binx1=binx-1
            binx2=binx
    x1=args.hbin.GetXaxis().GetBinCenter(binx1)
    x2=args.hbin.GetXaxis().GetBinCenter(binx2)

    ymin=args.hbin.GetYaxis().GetBinCenter(1)
    ymax=args.hbin.GetYaxis().GetBinCenter(args.hbin.GetYaxis().GetLast())
    if y<ymin:
        biny1=1
        biny2=2
    elif y>=ymax:
        biny1=args.hbin.GetYaxis().GetLast()-1
        biny2=args.hbin.GetYaxis().GetLast()
    else:
        biny=args.hbin.GetYaxis().FindBin(y)
        if y>=args.hbin.GetYaxis().GetBinCenter(biny):
            biny1=biny
            biny2=biny+1
        else:
            biny1=biny-1
            biny2=biny
    y1=args.hbin.GetYaxis().GetBinCenter(biny1)
    y2=args.hbin.GetYaxis().GetBinCenter(biny2)
    
    z11=args.zptweight[args.hbin.GetBin(binx1,biny1)].Eval(pt)
    z12=args.zptweight[args.hbin.GetBin(binx1,biny2)].Eval(pt)
    z21=args.zptweight[args.hbin.GetBin(binx2,biny1)].Eval(pt)
    z22=args.zptweight[args.hbin.GetBin(binx2,biny2)].Eval(pt)
    return ((x2-x)*(y2-y)*z11+(x2-x)*(y-y1)*z12+(x-x1)*(y2-y)*z12+(x-x1)*(y-y1)*z22)/(x2-x1)/(y2-y1)

def EvalZptWeight(args):
    npar=6
    mbins=args.hbin.GetXaxis().GetXbins()
    ybins=args.hbin.GetYaxis().GetXbins()
    for im in range(len(mbins)-1):
        for iy in range(len(ybins)-1):
            hdata=GetHist(args,0,mrange=[mbins[im],mbins[im+1]],yrange=[ybins[iy],ybins[iy+1]])
            hdy=GetHist(args,1,mrange=[mbins[im],mbins[im+1]],yrange=[ybins[iy],ybins[iy+1]],zptweight=True)
            hdata.Scale(1/hdata.Integral())
            hdy.Scale(1/hdy.Integral())
            hdata.Divide(hdy)
            for i in range(hdata.GetNcells()):
                sf=GetZptWeight(args,(mbins[im]+mbins[im+1])/2,(ybins[iy]+ybins[iy+1])/2,hdata.GetBinCenter(i))
                hdata.SetBinContent(i,sf*hdata.GetBinContent(i))
                hdata.SetBinError(i,sf*hdata.GetBinError(i))
            hdata.SetTitle("ratio_m{}_y{}".format(im,iy))
            hdata.Fit(args.zptweight[args.hbin.GetBin(im+1,iy+1)])
            hdata.Draw()
            hdata.GetYaxis().SetRangeUser(0.5,1.5)
    
    args.mPlotter.pdir.Delete()
    args.ePlotter.pdir.Delete()
    
def GetMinWidth(hist):
    minwidth=1e6
    for i in range(hist.GetNcells()):
        if minwidth>hist1.GetBinWidth(i): 
            minwidth=hist1.GetBinWidth(i)
    return minwidth

def ToFixedBinWidth(hist,width):
    nbins=int((hist.GetXaxis().GetXmax()-hist.GetXaxis().GetXmin())/width)
    out=rt.TH1D(hist.GetName()+"_fixed",hist.GetTitle()+"_fixed",nbins,hist.GetXaxis().GetXmin(),hist.GetXaxis().GetXmin()+nbins*width)
    for i in range(1,out.GetNcells()-1):
        high=out.GetBinLowEdge(i+1)
        low=out.GetBinLowEdge(i)
        hist_bin_start=hist.FindBin(low+width*1e-6)
        hist_bin_end=hist.FindBin(high-width*1e-6)
        val=0
        err2=0
        for j in range(hist_bin_start,hist_bin_end+1):
            this_val=hist.GetBinContent(j)
            this_err2=hist.GetBinError(j)*hist.GetBinError(j)*hist.GetBinWidth(j)/width
            weight=(min(hist.GetBinLowEdge(j+1),high)-max(hist.GetBinLowEdge(j),low))/width
            val+=weight*this_val
            err2+=weight*this_err2
        out.SetBinContent(i,val)
        out.SetBinError(i,sqrt(err2))
    return out

def Truncate(hist,prob=0.99,xmax=None):
    total=hist.Integral("width")
    integ=0
    bins=[]
    for i in range(1,hist.GetNcells()-1):
        bins+=[hist.GetBinLowEdge(i)]
        integ+=hist.GetBinContent(i)*hist.GetBinWidth(i)
        if integ/total>prob or (xmax is not None and xmax<=hist.GetBinLowEdge(i+1)): 
            bins+=[hist.GetBinLowEdge(i+1)]
            break
    out=rt.TH1D(hist.GetName()+"_truncated",hist.GetTitle()+"_truncated",len(bins)-1,array("d",bins))
    for i in range(1,hist.GetNcells()-1):
        out.SetBinContent(i,hist.GetBinContent(i))
        out.SetBinError(i,hist.GetBinError(i))
    return out

def SimpleConvolution(hist1,hist2,width=None):
    out=None
    if width is None: width=min(GetMinWidth(hist1),GetMinWidth(hist2))
    print hist1.Integral("width"),hist2.Integral("width")
    hist1=ToFixedBinWidth(hist1,width)
    hist2=ToFixedBinWidth(hist2,width)
    print hist1.Integral("width"),hist2.Integral("width")
    nbins=int((hist1.GetXaxis().GetXmax()+hist2.GetXaxis().GetXmax()-hist1.GetXaxis().GetXmin()-hist2.GetXaxis().GetXmin())/width)
    out=rt.TH1D(hist1.GetName()+"_"+hist2.GetName()+"_conv",hist1.GetTitle()+"_"+hist2.GetTitle()+"_conv",nbins,hist1.GetXaxis().GetXmin()+hist2.GetXaxis().GetXmin(),hist1.GetXaxis().GetXmin()+hist2.GetXaxis().GetXmin()+nbins*width)
    
    for i in range(1,hist1.GetNcells()-1):
        for j in range(1,hist2.GetNcells()-1):
            val=hist1.GetBinContent(i)*hist2.GetBinContent(j)*width
            err2=(hist1.GetBinContent(i)*hist2.GetBinError(j)*width)**2+(hist2.GetBinContent(i)*hist1.GetBinError(j)*width)**2
            ibin=out.FindBin(hist1.GetBinCenter(i)+hist2.GetBinCenter(j)-width*1e-6)
            out.SetBinContent(ibin,out.GetBinContent(ibin)+val)
            out.SetBinError(ibin,sqrt(out.GetBinError(ibin)**2+err2))
    return out

def test(args):
    SetupZptWeight(args)
    for i in range(10):
        print "iter",i
        EvalZptWeight(args)
    CheckZptWeight(args)
    
def run(args):
    cmd=""
    skim=" --skim SkimTree_Dilepton"
    cmdtemp="SKFlat.py -a ZptWeight -i {} -n {} -e {} --nmax 150 "
    runlist=[]
    for era in args.era:
        runlist+=[["WW_pythia"+skim,10,era]]
        runlist+=[["WZ_pythia"+skim,10,era]]
        runlist+=[["ZZ_pythia"+skim,10,era]]
        runlist+=[["WJets_MG"+skim,10,era]]
        runlist+=[["TTLL_powheg"+skim,30,era]]
        runlist+=[["SingleTop_tW_top_NoFullyHad"+skim,10,era]]
        runlist+=[["SingleTop_tW_antitop_NoFullyHad"+skim,10,era]]
        if args.dy=="MiNNLO":
            if "e" in ",".join(args.channel):
                runlist+=[["DYJetsToEE_MiNNLO",50,era]]
            if "m" in ",".join(args.channel):
                runlist+=[["DYJetsToMuMu_MiNNLO",50,era]]
            runlist+=[["DYJetsToTauTau_MiNNLO"+skim,10,era]]
        else:
            runlist+=[[args.dy,50,era]]
        if "e" in ",".join(args.channel):
            if era=="2018":
                runlist+=[["EGamma"+skim,20,era]]
            else:
                if "ee" in args.channel:
                    runlist+=[["DoubleEG"+skim,20,era]]
                if "el" in args.channel:
                    runlist+=[["SingleElectron"+skim,20,era]]
        if "mm" in args.channel:
            runlist+=[["DoubleMuon"+skim,20,era]]
        if "mu" in args.channel:
            runlist+=[["SingleMuon"+skim,20,era]]
            
    for a in runlist:
        cmd+=cmdtemp.format(a[0],a[1],a[2])+" & sleep 10;\n"
    cmd+="wait;\n"
    print(cmd)
    if args.dry==True: return;
    os.system(cmd)

if __name__=="__main__":
    

    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument("action",help="run, test")
    parser.add_argument("--dy",default="DYJets",help="DYJets, DYJets_MG, MiNNLO")
    parser.add_argument("--era",default="all",help="eras separated by commas. Available: 2016preVFP(2016a), 2016postVFP(2016b), 2017, 2018")
    parser.add_argument("--channel",default="mm,ee",help="channels separated by commas. Available: ee, el, mm, mu")
    parser.add_argument("--dry",default=False,action="store_true")
                       
    args=parser.parse_args()
    
    if args.era.lower() in ["all"]: 
        args.era="2016preVFP,2016postVFP,2017,2018"
        #args.era="2017"
        if args.dy=="MiNNLO": args.era="2016preVFP,2016postVFP"
    args.era=args.era.replace("2016a","2016preVFP").replace("2016b","2016postVFP").split(",")
    args.channel=args.channel.split(",")

    if args.dy=="DYJets":
        args.mPlotter=rt.AFBPlotter("data-tau_amc-vv-wjets-tttw-1.7*ss_amc amc","ZptWeight")
        args.ePlotter=rt.AFBPlotter("data-tau_amc-vv-wjets-tttw-ss_amc amc","ZptWeight")
    elif args.dy=="DYJets_MG":
        args.mPlotter=rt.AFBPlotter("data-tau_mg-vv-wjets-tttw-1.7*ss_mg mg","ZptWeight")
        args.ePlotter=rt.AFBPlotter("data-tau_mg-vv-wjets-tttw-ss_mg mg","ZptWeight")
    elif args.dy=="MiNNLO":
        args.mPlotter=rt.AFBPlotter("data-tau_mi-vv-wjets-tttw-1.7*ss_mi mi","ZptWeight")
        args.ePlotter=rt.AFBPlotter("data-tau_mi-vv-wjets-tttw-ss_mi mi","ZptWeight")
    rt.Verbosity=0

    if args.action in ["run"]:
        run(args)
    elif args.action in ["test"]:
        test(args)
    else:
        print "unavailable action",args.action
        exit(1)
    
