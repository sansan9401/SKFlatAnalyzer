import os
import scipy.integrate as integrate
from scipy.optimize import minimize
from math import exp,sqrt,log
from array import array
import ROOT as rt
rt.gROOT.SetBatch(True)
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

    name="pt"
    if ientry==0: name="data_"+name
    else: name="dy_"+name
    if gen: name="gen_"+name
    if zptweight: name+="_zpt"
    hist1d=hist.ProjectionZ(name,ixmin,ixmax,iymin,iymax)
    hist1d.SetDirectory(0)
    hist.Delete()
    return hist1d

def FuncZptWeightG(xx,par):
    #npar_=14
    kk=[0,2,4,7,10,15,20,30,40,60,100,200,400]
    x=xx[0]
    aa=[par[i] for i in range(len(kk))]
    bb=[None]*(len(kk)-1)+[par[len(kk)]]
    cc=[None]*(len(kk)-1)+[0]
    dd=[None]*(len(kk)-1)+[0]
    for i in range(len(kk)-1,-1,-1):
        if x>=kk[i] or i==0:
            return aa[i]+bb[i]*(x-kk[i])+cc[i]*(x-kk[i])**2+dd[i]*(x-kk[i])**3
        r=kk[i]-kk[i-1]
        dd[i-1]=cc[i]/r-bb[i]/r**2+(aa[i]-aa[i-1])/r**3
        cc[i-1]=cc[i]-3*dd[i-1]*r
        bb[i-1]=bb[i]-2*cc[i-1]*r-3*dd[i-1]*r**2
    return None

def FuncZptWeightY(xx,par):
    #npar_y=7
    kk=[0,10,20,40,60,200]
    x=xx[0]
    aa=[par[i] for i in range(len(kk))]
    bb=[None]*(len(kk)-1)+[par[len(kk)]]
    cc=[None]*(len(kk)-1)+[0]
    dd=[None]*(len(kk)-1)+[0]
    for i in range(len(kk)-1,-1,-1):
        if x>=kk[i] or i==0:
            return aa[i]+bb[i]*(x-kk[i])+cc[i]*(x-kk[i])**2+dd[i]*(x-kk[i])**3
        r=kk[i]-kk[i-1]
        dd[i-1]=cc[i]/r-bb[i]/r**2+(aa[i]-aa[i-1])/r**3
        cc[i-1]=cc[i]-3*dd[i-1]*r
        bb[i-1]=bb[i]-2*cc[i-1]*r-3*dd[i-1]*r**2
    return None

def FuncZptWeightM(xx,par):
    #npar_m=8
    kk=[0,5,10,20,30,55,100]
    x=xx[0]
    aa=[par[i] for i in range(len(kk))]
    bb=[None]*(len(kk)-1)+[par[len(kk)]]
    cc=[None]*(len(kk)-1)+[0]
    dd=[None]*(len(kk)-1)+[0]
    for i in range(len(kk)-1,-1,-1):
        if x>=kk[i] or i==0:
            return aa[i]+bb[i]*(x-kk[i])+cc[i]*(x-kk[i])**2+dd[i]*(x-kk[i])**3
        r=kk[i]-kk[i-1]
        dd[i-1]=cc[i]/r-bb[i]/r**2+(aa[i]-aa[i-1])/r**3
        cc[i-1]=cc[i]-3*dd[i-1]*r
        bb[i-1]=bb[i]-2*cc[i-1]*r-3*dd[i-1]*r**2
    return None

def Scale(func,sf):
    if not func.InheritsFrom("TF1"):
        print "cannot scaling "+func.ClassName()
        return
    for i in range(func.GetNpar()):
        func.SetParameter(i,func.GetParameter(i)*sf)
        func.SetParError(i,func.GetParError(i)*sf)

def SetupZptWeight(args):
    args.npar_g=14
    args.npar_y=7
    args.npar_m=8

    if args.input:
        args.input=rt.TFile(args.input)
            
    ## step 1 global correction
    args.zptweight_g=rt.TF1("zptweight_g",FuncZptWeightG,0,650,args.npar_g)
    args.zptweight_g.SetNpx(1000)
    if args.input:
        args.zptweight_g.SetParameters(args.input.Get("zptweight_g").GetParameters())
    else:
        for i in range(args.npar_g-1):
            args.zptweight_g.SetParameter(i,1)

    ## step 2 rapidity-dependent correction
    if args.input:
        args.yaxis=args.input.Get("yaxis")
    else:
        ybins=[0.0,0.4,1.0,1.6,2.4]
        args.yaxis=rt.TAxis(len(ybins)-1,array("d",ybins))
        args.yaxis.SetName("yaxis")

    args.zptweight_y=[None]*(args.yaxis.GetNbins()+2)
    for i in range(1,args.yaxis.GetNbins()+1):
        args.zptweight_y[i]=rt.TF1("zptweight_y{}".format(i),FuncZptWeightY,0,650,args.npar_y)
        args.zptweight_y[i].SetNpx(1000)
        if args.input:
            args.zptweight_y[i].SetParameters(args.input.Get("zptweight_y{}".format(i)).GetParameters())
        else:
            for ip in range(args.npar_y-1):
                args.zptweight_y[i].SetParameter(ip,1)
        args.zptweight_y[i].FixParameter(args.npar_y-1,0)

    ## step 3 mass-dependent correction
    #args.input=None ## FIXME temp for m study
    if args.input:
        args.maxis=args.input.Get("maxis")
    else:
        mbins=[52,80,90,100,150]
        args.maxis=rt.TAxis(len(mbins)-1,array("d",mbins))
        args.maxis.SetName("maxis")
    args.zptweight_m=[None]*(args.maxis.GetNbins()+2)
    for i in range(1,args.maxis.GetNbins()+1):
        args.zptweight_m[i]=rt.TF1("zptweight_m{}".format(i),FuncZptWeightM,0,650,args.npar_m)
        args.zptweight_m[i].SetNpx(1000)
        if args.input: 
            args.zptweight_m[i].SetParameters(args.input.Get("zptweight_m{}".format(i)).GetParameters())
        else:
            for ip in range(args.npar_m-1):
                args.zptweight_m[i].SetParameter(ip,1)
        args.zptweight_m[i].FixParameter(args.npar_m-1,0)

def CheckZptWeight(args):
    ms=[52,80,100,150]
    ys=[0,0.8,1.6,2.4]
    for im in range(len(ms)):
        c=rt.TCanvas()
        gs=[]
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
        args.out.cd()
        c.Write("check_zpt_m{}".format(ms[im]))

    for iy in range(len(ys)):
        c=rt.TCanvas()
        gs=[]
        for im in range(len(ms)):
            npoints=400
            xpoints=[i for i in range(npoints)]
            ypoints=[GetZptWeight(args,ms[im],ys[iy],x) for x in xpoints]
            gs+=[rt.TGraph(npoints,array("d",xpoints),array("d",ypoints))]
            if im==0:
                gs[-1].Draw("al")
                gs[-1].SetTitle("zptweight at y={}".format(ys[iy]))
                gs[-1].GetYaxis().SetRangeUser(0.5,1.5)
            else:
                gs[-1].SetLineColor(im+1)
                gs[-1].Draw("same l")
        args.out.cd()
        c.Write("check_zpt_y{}".format(ys[iy]))

def CheckResidue(args,prefix=""):
    ## global
    hdata=GetHist(args,0,mrange=[80,100],yrange=[0,2.4])
    hdy=GetHist(args,1,mrange=[80,100],yrange=[0,2.4],zptweight=True)
    hdata.Scale(1/hdata.Integral())
    hdy.Scale(1/hdy.Integral())
    hdata.Divide(hdy)
    hdata.GetYaxis().SetRangeUser(0.5,1.5)
    args.out.cd()
    hdata.Write(prefix+"residue_g")

    ## rapidity
    yaxis=args.yaxis
    for iy in range(1,yaxis.GetNbins()+1):
        hdata=GetHist(args,0,mrange=[80,100],yrange=[yaxis.GetBinLowEdge(iy),yaxis.GetBinUpEdge(iy)])
        hdy=GetHist(args,1,mrange=[80,100],yrange=[yaxis.GetBinLowEdge(iy),yaxis.GetBinUpEdge(iy)],zptweight=True)
        hdata.Scale(1/hdata.Integral())
        hdy.Scale(1/hdy.Integral())
        hdata.Divide(hdy)
        hdata.GetYaxis().SetRangeUser(0.5,1.5)
        args.out.cd()
        hdata.Write(prefix+"residue_y{}".format(iy))
    
    ## mass
    maxis=args.maxis
    for im in range(1,maxis.GetNbins()+1):
        hdata=GetHist(args,0,mrange=[maxis.GetBinLowEdge(im),maxis.GetBinUpEdge(im)])
        hdy=GetHist(args,1,mrange=[maxis.GetBinLowEdge(im),maxis.GetBinUpEdge(im)],zptweight=True)
        hdata.Scale(1/hdata.Integral())
        hdy.Scale(1/hdy.Integral())
        hdata.Divide(hdy)
        hdata.GetYaxis().SetRangeUser(0.5,1.5)
        args.out.cd()
        hdata.Write(prefix+"residue_m{}".format(im))
    

def GetZptWeight(args,mass,rapidity,pt,step="GYM"):
    x=mass
    y=abs(rapidity)
    sf=1.
    step=step.upper()

    ## step 1 (G) global correction
    if "G" in step:
        sf*=args.zptweight_g.Eval(pt)

    ## step 2 (Y) rapidity-depenent correction
    if "Y" in step:
        ymin=args.yaxis.GetBinCenter(1)
        ymax=args.yaxis.GetBinCenter(args.yaxis.GetNbins())
        if y<ymin:
            biny1=1
            biny2=2
        elif y>=ymax:
            biny1=args.yaxis.GetNbins()-1
            biny2=args.yaxis.GetNbins()
        else:
            biny=args.yaxis.FindBin(y)
            if y>=args.yaxis.GetBinCenter(biny):
                biny1=biny
                biny2=biny+1
            else:
                biny1=biny-1
                biny2=biny
        y1=args.yaxis.GetBinCenter(biny1)
        y2=args.yaxis.GetBinCenter(biny2)
        sf*=( (y2-y)*args.zptweight_y[biny1].Eval(pt) + (y-y1)*args.zptweight_y[biny2].Eval(pt) )/(y2-y1)
    
    ## step 3 mass-depenent correction
    if "M" in step:
        xmin=args.maxis.GetBinCenter(1)
        xmax=args.maxis.GetBinCenter(args.maxis.GetNbins())
        if x<xmin:
            binx1=1
            binx2=2
        elif x>=xmax:
            binx1=args.maxis.GetNbins()-1
            binx2=args.maxis.GetNbins()
        else:
            binx=args.maxis.FindBin(x)
            if x>=args.maxis.GetBinCenter(binx):
                binx1=binx
                binx2=binx+1
            else:
                binx1=binx-1
                binx2=binx
        x1=args.maxis.GetBinCenter(binx1)
        x2=args.maxis.GetBinCenter(binx2)
        sf*=( (x2-x)*args.zptweight_m[binx1].Eval(pt) + (x-x1)*args.zptweight_m[binx2].Eval(pt) )/(x2-x1)
    return sf

def EvalZptWeight(args):
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

def EvalZptWeightG(args,iteration=1):
    for it in range(iteration):
        print "[EvalZptWeightG] iter {}".format(it)
        hdata=GetHist(args,0,mrange=[80,100],yrange=[0,2.4])
        hdy=GetHist(args,1,mrange=[80,100],yrange=[0,2.4],zptweight=True)
        hdata.Scale(1/hdata.Integral())
        hdy.Scale(1/hdy.Integral())
        hdata.Divide(hdy)
        hdata.GetYaxis().SetRangeUser(0.5,1.5)
        args.out.cd()
        hdata.Write("gresidue_iter{}".format(it))
        for i in range(hdata.GetNcells()):
            sf=GetZptWeight(args,91,0,hdata.GetBinCenter(i))
            hdata.SetBinContent(i,sf*hdata.GetBinContent(i))
            hdata.SetBinError(i,sf*hdata.GetBinError(i))
        c=rt.TCanvas()
        hdata.SetTitle("gfit_iter{}".format(it))
        hdata.Fit(args.zptweight_g)
        hdata.Draw()
        args.out.cd()
        c.Write("gfit_iter{}".format(it))
        args.mPlotter.pdir.Delete()
        args.ePlotter.pdir.Delete()
    
    hdy_nozptweight=GetHist(args,1,zptweight=False,gen=True)
    hdy=GetHist(args,1,zptweight=True,gen=True)
    normsf=hdy_nozptweight.Integral()/hdy.Integral()
    print "[EvalZptWeightG] norm sf={}".format(normsf)
    Scale(args.zptweight_g,normsf)
    hdata=GetHist(args,0,mrange=[80,100],yrange=[0,2.4])
    hdy=GetHist(args,1,mrange=[80,100],yrange=[0,2.4],zptweight=True)
    hdata.Scale(1/hdata.Integral())
    hdy.Scale(1/hdy.Integral())
    hdata.Divide(hdy)
    hdata.GetYaxis().SetRangeUser(0.5,1.5)
    args.out.cd()
    hdata.Write("gresidue_norm")

def EvalZptWeightY(args,iteration=1):
    yaxis=args.yaxis
    for it in range(iteration):
        for iy in range(1,yaxis.GetNbins()+1):
            print "[EvalZptWeightY] iter {} y {}".format(it,iy)
            hdata=GetHist(args,0,mrange=[80,100],yrange=[yaxis.GetBinLowEdge(iy),yaxis.GetBinUpEdge(iy)])
            hdy=GetHist(args,1,mrange=[80,100],yrange=[yaxis.GetBinLowEdge(iy),yaxis.GetBinUpEdge(iy)],zptweight=True)
            hdata.Scale(1/hdata.Integral())
            hdy.Scale(1/hdy.Integral())
            hdata.Divide(hdy)
            hdata.GetYaxis().SetRangeUser(0.5,1.5)
            args.out.cd()
            hdata.Write("yresidue_iter{}_y{}".format(it,iy))
            for i in range(hdata.GetNcells()):
                sf=GetZptWeight(args,91,yaxis.GetBinCenter(iy),hdata.GetBinCenter(i),step="Y")
                hdata.SetBinContent(i,sf*hdata.GetBinContent(i))
                hdata.SetBinError(i,sf*hdata.GetBinError(i))
            c=rt.TCanvas()
            hdata.SetTitle("yfit_iter{}_y{}".format(it,iy))
            hdata.Fit(args.zptweight_y[iy])
            hdata.Draw()
            args.out.cd()
            c.Write("yfit_iter{}_y{}".format(it,iy))
            args.mPlotter.pdir.Delete()
            args.ePlotter.pdir.Delete()
        
    for iy in range(1,yaxis.GetNbins()+1):
        hdy_nozptweight=GetHist(args,1,yrange=[yaxis.GetBinLowEdge(iy),yaxis.GetBinUpEdge(iy)],zptweight=False,gen=True)
        hdy=GetHist(args,1,yrange=[yaxis.GetBinLowEdge(iy),yaxis.GetBinUpEdge(iy)],zptweight=True,gen=True)
        normsf=hdy_nozptweight.Integral()/hdy.Integral()
        print "[EvalZptWeightY] ybin {} norm sf={}".format(iy,normsf)
        Scale(args.zptweight_y[iy],normsf)

def EvalZptWeightM(args,iteration=1):
    maxis=args.maxis
    for it in range(iteration):
        for im in range(1,maxis.GetNbins()+1):
            print "[EvalZptWeightM] iter {} m {}".format(it,im)
            hdata=GetHist(args,0,mrange=[maxis.GetBinLowEdge(im),maxis.GetBinUpEdge(im)],yrange=[0,2.4])
            hdy=GetHist(args,1,mrange=[maxis.GetBinLowEdge(im),maxis.GetBinUpEdge(im)],yrange=[0,2.4],zptweight=True)
            hdata.Scale(1/hdata.Integral())
            hdy.Scale(1/hdy.Integral())
            hdata.Divide(hdy)
            hdata.GetYaxis().SetRangeUser(0.5,1.5)
            args.out.cd()
            hdata.Write("mresidue_iter{}_m{}".format(it,im))
            for i in range(hdata.GetNcells()):
                sf=GetZptWeight(args,maxis.GetBinCenter(im),0,hdata.GetBinCenter(i),step="M")
                hdata.SetBinContent(i,sf*hdata.GetBinContent(i))
                hdata.SetBinError(i,sf*hdata.GetBinError(i))
            c=rt.TCanvas()
            hdata.SetTitle("mfit_iter{}_m{}".format(it,im))
            args.zptweight_m[im].FixParameter(0,1)
            args.zptweight_m[im].FixParameter(args.npar_m-1,0)
            hdata.Fit(args.zptweight_m[im])
            hdata.Draw()
            args.out.cd()
            c.Write("mfit_iter{}_m{}".format(it,im))
            args.mPlotter.pdir.Delete()
            args.ePlotter.pdir.Delete()
        
    for im in range(1,maxis.GetNbins()+1):
        hdy_nozptweight=GetHist(args,1,mrange=[maxis.GetBinLowEdge(im),maxis.GetBinUpEdge(im)],zptweight=False,gen=True)
        hdy=GetHist(args,1,mrange=[maxis.GetBinLowEdge(im),maxis.GetBinUpEdge(im)],zptweight=True,gen=True)
        normsf=hdy_nozptweight.Integral()/hdy.Integral()
        print "[EvalZptWeightM] mbin {} norm sf={}".format(im,normsf)
        Scale(args.zptweight_m[im],normsf)
    
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

def SaveZptWeight(args):
    args.out.cd()
    args.yaxis.Write()
    args.maxis.Write()
    args.zptweight_g.Write()
    for f in args.zptweight_y:
        if f: f.Write()
    for f in args.zptweight_m:
        if f: f.Write()

def test(args):
    SetupZptWeight(args)
    EvalZptWeightG(args,iteration=5)
    EvalZptWeightY(args,iteration=5)
    EvalZptWeightM(args,iteration=5)
    CheckResidue(args,prefix="after_")
    CheckZptWeight(args)
    SaveZptWeight(args)
    
def run(args):
    cmd=""
    skim=" --skim SkimTree_Dilepton"
    cmdtemp="SKFlat.py -a ZptWeight -i {} -n {} -e {} --nmax 150 "
    runlist=[]
    for era in args.era:
        if "bg" in args.run:
            runlist+=[["WW_pythia"+skim,10,era]]
            runlist+=[["WZ_pythia"+skim,10,era]]
            runlist+=[["ZZ_pythia"+skim,10,era]]
            runlist+=[["WJets_MG"+skim,10,era]]
            runlist+=[["TTLL_powheg"+skim,30,era]]
            runlist+=[["SingleTop_tW_top_NoFullyHad"+skim,10,era]]
            runlist+=[["SingleTop_tW_antitop_NoFullyHad"+skim,10,era]]
        if "dy" in args.run:
            if args.dy=="MiNNLO":
                if "e" in ",".join(args.channel):
                    runlist+=[["DYJetsToEE_MiNNLO",50,era]]
                if "m" in ",".join(args.channel):
                    runlist+=[["DYJetsToMuMu_MiNNLO",50,era]]
                runlist+=[["DYJetsToTauTau_MiNNLO"+skim,10,era]]
            else:
                runlist+=[[args.dy,50,era]]
        if "data" in args.run:
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
    parser.add_argument("--run",default="data,dy,bg",help="run list separated by commas. Available: data, dy,bg")
    parser.add_argument("--in",dest="input",default=None,help="input file path")
    parser.add_argument("--out",default="zptout.root",help="out file path")
    parser.add_argument("--dry",default=False,action="store_true")
                       
    args=parser.parse_args()
    
    if args.era.lower() in ["all"]: 
        args.era="2016preVFP,2016postVFP,2017,2018"
        #args.era="2017"
        if args.dy=="MiNNLO": args.era="2016preVFP,2016postVFP"
    args.era=args.era.replace("2016a","2016preVFP").replace("2016b","2016postVFP").split(",")
    args.channel=args.channel.split(",")
    args.run=args.run.split(",")

    if args.out==args.input:
        print "cannot use out==in"
        exit(1)
    args.out=rt.TFile(args.out,"recreate")

    rt.Verbosity=1
    if args.dy=="DYJets":
        args.mPlotter=rt.AFBPlotter("data-tau_amc-vv-wjets-tttw-1.7*ss_amc amc","ZptWeight")
        args.ePlotter=rt.AFBPlotter("data-tau_amc-vv-wjets-tttw-ss_amc amc","ZptWeight")
    elif args.dy=="DYJets_MG":
        args.mPlotter=rt.AFBPlotter("data-tau_mg-vv-wjets-tttw-1.7*ss_mg mg","ZptWeight")
        args.ePlotter=rt.AFBPlotter("data-tau_mg-vv-wjets-tttw-ss_mg mg","ZptWeight")
    elif args.dy=="MiNNLO":
        args.mPlotter=rt.AFBPlotter("data-tau_mi-vv-wjets-tttw-1.7*ss_mi mi","ZptWeight")
        args.ePlotter=rt.AFBPlotter("data-tau_mi-vv-wjets-tttw-ss_mi mi","ZptWeight")

    if args.action in ["run"]:
        run(args)
    elif args.action in ["test"]:
        test(args)
    else:
        print "unavailable action",args.action
        exit(1)
    
