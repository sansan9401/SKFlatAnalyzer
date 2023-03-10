import os,ctypes,array,copy
import ROOT
#ROOT.gROOT.SetBatch(True)
ROOT.gROOT.LoadMacro("./Plotter/SKFlatPlotter.cc")
ROOT.TH1.AddDirectory(0)

def GetHist(filename,histname):
    f=ROOT.TFile(filename)
    h=f.Get(histname)
    if not h: print "no {} in {}".format(histname,filename)
    return h

def GetEra(era):
    if era=="2016a": 
        return "2016preVFP"
    elif era=="2016b": 
        return "2016postVFP"
    else: 
        return era

def GetEraShort(era):
    if era=="2016preVFP": 
        return "2016a"
    elif era=="2016postVFP": 
        return "2016b"
    else: 
        return era
    

def Rebin2D(hist,xbins,ybins):
    histnew=ROOT.TH2D(hist.GetName(),hist.GetTitle(),len(xbins)-1,array.array("d",xbins),len(ybins)-1,array.array("d",ybins))
    histnew.SetStats(0)
    for i in range(hist.GetNcells()):
        val=hist.GetBinContent(i)
        err=hist.GetBinError(i)
        binx=ctypes.c_int()
        biny=ctypes.c_int()
        binz=ctypes.c_int()
        hist.GetBinXYZ(i,binx,biny,binz)
        binx=binx.value
        biny=biny.value
        x=hist.GetXaxis().GetBinCenter(binx)
        y=hist.GetYaxis().GetBinCenter(biny)
        
        j=histnew.FindBin(x,y)
        histnew.SetBinContent(j,histnew.GetBinContent(j)+val)
        histnew.SetBinError(j,(histnew.GetBinError(j)**2+err**2)**0.5)
    return histnew


def MakeTF(args,num,den):
    if not num: return None
    if not den: return None
    hnum=Rebin2D(num,args.tfxbins,args.tfybins)
    hden=Rebin2D(den,args.tfxbins,args.tfybins)
    hnum.Divide(hden)
    return hnum

def MakeSet(hists,name):
    print "[MakeSet]", name
    rt=[]
    if not hists[0][0]: return rt
    hist=hists[0][0].Clone(name)
    hist.SetTitle(name)
    for i in range(1,len(hists)):
        for ibin in range(hist.GetNcells()):
            val=hist.GetBinContent(ibin)
            err=hist.GetBinError(ibin)
            maxdiff=0
            for j in range(len(hists[i])):
                diff=abs(hists[i][j].GetBinContent(ibin)-val)
                if diff>maxdiff:
                    maxdiff=diff
            hist.SetBinError(ibin,(err**2+maxdiff**2)**0.5)
    rt+=[hist]
    for i in range(len(hists)):
        for j in range(len(hists[i])):
            hists[i][j].SetName(name+"_s{}m{}".format(i,j))
            hists[i][j].SetTitle(name+"_s{}m{}".format(i,j))
            rt+=[hists[i][j].Clone()]
    return rt

if __name__=="__main__":
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument("-e","--era",default="2016a,2016b,2017,2018",type=str)
    parser.add_argument("-c","--channel",default="ee,mm,ej,mj")
    parser.add_argument("-r","--region",default=",noZ/,0bjet/noZ/,nbjet/noZ/,cpt/noZ/,cpt/0bjet/noZ/,cpt/nbjet/noZ/,mpt/noZ/,mpt/0bjet/noZ/,mpt/nbjet/noZ/")
    args_raw=parser.parse_args()

    args_raw.era=args_raw.era.split(",")
    args_raw.channel=args_raw.channel.split(",")
    args_raw.region=args_raw.region.split(",")

    save=[]
    for era in args_raw.era:
        for channel in args_raw.channel:
            for region in args_raw.region:
                args=copy.deepcopy(args_raw)
                args.plotter=ROOT.SKFlatPlotter("data-mi-tau_mi-wjets-vv-tt-st-aa","FakeAnalyzer")
                args.era=GetEra(era)
                args.erashort=GetEraShort(era)
                if args.era not in ["2016preVFP","2016postVFP","2017","2018"]:
                    print "Unknown era {}".format(args.era)
                    continue
                args.region=region
                args.channel=channel
                args.tfxbins=[0.,1.5,2.5]
                args.tfybins=[10,12.5,15,17.5,20,22.5,25,27.5,30,35,40,50,100,1000]
                if args.channel=="ee":
                    args.l0tf_num_str="eE{}/{}/ss_l0etapt".format(args.erashort,region)
                    args.l0tf_den_str="EE{}/{}/ss_l0etapt".format(args.erashort,region)
                    args.l1tf_num_str="Ee{}/{}/ss_l1etapt".format(args.erashort,region)
                    args.l1tf_den_str="EE{}/{}/ss_l1etapt".format(args.erashort,region)

                elif args.channel=="mm":
                    args.l0tf_num_str="mM{}/{}/ss_l0etapt".format(args.erashort,region)
                    args.l0tf_den_str="MM{}/{}/ss_l0etapt".format(args.erashort,region)
                    args.l1tf_num_str="Mm{}/{}/ss_l1etapt".format(args.erashort,region)
                    args.l1tf_den_str="MM{}/{}/ss_l1etapt".format(args.erashort,region)

                elif args.channel=="ej":
                    args.l0tf_num_str="ej{}/{}/l0etapt".format(args.erashort,region)
                    args.l0tf_den_str="Ej{}/{}/l0etapt".format(args.erashort,region)
                    args.l1tf_num_str="ej{}/{}/l0etapt".format(args.erashort,region)
                    args.l1tf_den_str="Ej{}/{}/l0etapt".format(args.erashort,region)
                    
                elif args.channel=="mj":
                    args.l0tf_num_str="mj{}/{}/l0etapt".format(args.erashort,region)
                    args.l0tf_den_str="Mj{}/{}/l0etapt".format(args.erashort,region)
                    args.l1tf_num_str="mj{}/{}/l0etapt".format(args.erashort,region)
                    args.l1tf_den_str="Mj{}/{}/l0etapt".format(args.erashort,region)

                else:
                    print "Unknown channel {}".format(args.channel)
                    continue

                hl0tf_s0m0=MakeTF(args,
                                  args.plotter.GetHist(0,args.l0tf_num_str,"noproject"),
                                  args.plotter.GetHist(0,args.l0tf_den_str,"noproject"))
                hl0tf_s1m0=MakeTF(args,
                                  args.plotter.GetHist(0,args.l0tf_num_str,"noproject scale:1.1:data"),
                                  args.plotter.GetHist(0,args.l0tf_den_str,"noproject scale:1.1:data"))
                hl0tf_s1m1=MakeTF(args,
                                  args.plotter.GetHist(0,args.l0tf_num_str,"noproject scale:0.9:data"),
                                  args.plotter.GetHist(0,args.l0tf_den_str,"noproject scale:0.9:data"))
                save+=MakeSet([[hl0tf_s0m0],[hl0tf_s1m0,hl0tf_s1m1]],"{}{}_{}l0tf".format(args.channel,args.era,args.region.replace("/","_")))

                hl1tf_s0m0=MakeTF(args,
                                  args.plotter.GetHist(0,args.l1tf_num_str,"noproject"),
                                  args.plotter.GetHist(0,args.l1tf_den_str,"noproject"))
                hl1tf_s1m0=MakeTF(args,
                                  args.plotter.GetHist(0,args.l1tf_num_str,"noproject scale:1.1:data"),
                                  args.plotter.GetHist(0,args.l1tf_den_str,"noproject scale:1.1:data"))
                hl1tf_s1m1=MakeTF(args,
                                  args.plotter.GetHist(0,args.l1tf_num_str,"noproject scale:0.9:data"),
                                  args.plotter.GetHist(0,args.l1tf_den_str,"noproject scale:0.9:data"))
                save+=MakeSet([[hl1tf_s0m0],[hl1tf_s1m0,hl1tf_s1m1]],"{}{}_{}l1tf".format(args.channel,args.era,args.region.replace("/","_")))

                        
    f=ROOT.TFile("FakeTF.root","recreate")
    for h in save:
        for i in range(h.GetNcells()):
            val=h.GetBinContent(i)
            if val<0:
                h.SetBinContent(i,0.)
                h.SetBinError(i,0.)
        h.Write()
    exit(0)
