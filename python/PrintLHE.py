import os,sys
import argparse
import ROOT

parser = argparse.ArgumentParser(description='SKFlat Command')
parser.add_argument('-i', dest='InputSample', default="")
parser.add_argument('-e', dest='Era', default="2017",help="2016preVFP(2016a), 2016postVFP(2016b), 2017, 2018")
parser.add_argument('-n', dest='n', default=10, type=int)
parser.add_argument('--skim', dest='Skim', default="", help="ex) SkimTree_Dilepton")
args = parser.parse_args()

if args.Era=="2016a": args.Era="2016preVFP"
if args.Era=="2016b": args.Era="2016postVFP"

SkimString = ""
if args.Skim!="":
  SkimString = args.Skim+"_"

SKFlat_WD = os.environ['SKFlat_WD']
SKFlatV = os.environ['SKFlatV']
SAMPLE_DATA_DIR = SKFlat_WD+'/data/'+SKFlatV+'/'+args.Era+'/Sample/'
filelistpath = SAMPLE_DATA_DIR+'/ForSNU/'+SkimString+args.InputSample+'.txt'
if not os.path.exists(filelistpath):
    print "No file "+filelistpath
    exit(1)
filepaths = os.popen("sed 's/#.*//' "+filelistpath+"|grep '.root'").read().split()
analyzer=ROOT.AnalyzerCore()
analyzer.SetTreeName("recoTree/SKFlat")
analyzer.MCSample=args.InputSample
analyzer.SetEra(args.Era)
analyzer.AddFile(filepaths[0])
analyzer.Init()
analyzer.initializeAnalyzer()
analyzer.initializeAnalyzerTools()
print "========= file:"+filepaths[0]+" ============="
for i in range(0,args.n):
    print "----------- entry "+str(i)+" -------------"
    analyzer.GetEntry(i)
    analyzer.beginEvent()
    lhes=analyzer.GetLHEs()
    for lhe in lhes:
        lhe.Print()
    analyzer.endEvent()

