#!/usr/bin/env python

import os,sys,time
import argparse
import datetime
from CheckJobStatus import *
from GetXSECTable import *
from TimeTools import *
import random
import subprocess

## Arguments

parser = argparse.ArgumentParser(description='SKFlat Command')
parser.add_argument('-a', dest='Analyzer', default="")
parser.add_argument('-i', dest='InputSample', default="")
parser.add_argument('-p', dest='DataPeriod', default="ALL")
parser.add_argument('-l', dest='InputSampleList', default="")
parser.add_argument('-n', dest='NJobs', default=1, type=int)
parser.add_argument('-o', dest='Outputdir', default="")
parser.add_argument('-q', dest='Queue', default="fastq")
parser.add_argument('-e', dest='Era', default="2017",help="2016preVFP(2016a), 2016postVFP(2016b), 2017, 2018")
parser.add_argument('-y', dest='Year', default="",help="deprecated. use -e")
parser.add_argument('--skim', dest='Skim', default="", help="ex) SkimTree_Dilepton")
parser.add_argument('--no_exec', action='store_true')
parser.add_argument('--FastSim', action='store_true')
parser.add_argument('--userflags', dest='Userflags', default="")
parser.add_argument('--nmax', dest='NMax', default=0, type=int, help="maximum running jobs")
parser.add_argument('--reduction', dest='Reduction', default=1, type=float)
parser.add_argument('--memory', dest='Memory', default=0, type=float)
parser.add_argument('--batchname',dest='BatchName', default="")
args = parser.parse_args()

if args.Year!="":
  print("-y is deprecated. Using -e (Era) instead")
  args.Era=args.Year

if args.Era=="2016a": args.Era="2016preVFP"
if args.Era=="2016b": args.Era="2016postVFP"

## make userflags as a list
Userflags = []
if args.Userflags != "":
  Userflags = (args.Userflags).split(',')

## Add Abosolute path for outputdir
if args.Outputdir!='':
  if args.Outputdir[0]!='/':
    args.Outputdir = os.getcwd()+'/'+args.Outputdir

## TimeStamp

# 1) dir/file name style
JobStartTime = datetime.datetime.now()
timestamp =  JobStartTime.strftime('%Y_%m_%d_%H%M%S')
# 2) log style
JobStartTime = datetime.datetime.now()
string_JobStartTime =  JobStartTime.strftime('%Y-%m-%d %H:%M:%S')
string_ThisTime = ""

## Environment Variables

USER = os.environ['USER']
if os.path.exists('python/UserInfo_'+USER+'.py'):
  exec('from UserInfo_'+USER+' import *')
else:
  print("No UserInfo file")
  exit(1)
SKFlatLogEmail = UserInfo['SKFlatLogEmail']
SKFlatLogWebDir = UserInfo['SKFlatLogWebDir']
LogEvery = UserInfo['LogEvery']

SCRAM_ARCH = os.environ['SCRAM_ARCH']
cmsswrel = os.environ['cmsswrel']
SKFlat_WD = os.environ['SKFlat_WD']
SKFlatV = os.environ['SKFlatV']
SAMPLE_DATA_DIR = SKFlat_WD+'/data/'+SKFlatV+'/'+args.Era+'/Sample/'
SKFlatRunlogDir = os.environ['SKFlatRunlogDir']
SKFlatOutputDir = os.environ['SKFlatOutputDir']
SKFlat_LIB_PATH = os.environ['SKFlat_LIB_PATH']
UID = str(os.getuid())
HOSTNAME = os.environ['HOSTNAME']
SampleHOSTNAME = HOSTNAME

## Check joblog email

if SKFlatLogEmail=='':
  print '[SKFlat.py] Put your email address in setup.sh'
  exit()
SendLogToWeb = True
if SKFlatLogWebDir=='':
  SendLogToWeb = False

## Check hostname

IsKISTI = ("sdfarm.kr" in HOSTNAME)
IsUI10 = ("ui10.sdfarm.kr" in HOSTNAME)
IsUI20 = ("ui20.sdfarm.kr" in HOSTNAME)
IsTAMSA1 = ("tamsa1" in HOSTNAME)
IsTAMSA2 = ("tamsa2" in HOSTNAME)
IsTAMSA = IsTAMSA1 or IsTAMSA2
IsKNU = ("knu" in HOSTNAME)
if IsKISTI:
  HOSTNAME = "KISTI"
  SampleHOSTNAME = "KISTI"
if IsTAMSA:
  if IsTAMSA1:
    HOSTNAME = "TAMSA1"
  elif IsTAMSA2:
    HOSTNAME = "TAMSA2"
  SampleHOSTNAME = "SNU"
if IsKNU:
  HOSTNAME = "KNU"
  SampleHOSTNAME = "KNU"

## Are you skimming trees?

IsSkimTree = "SkimTree" in args.Analyzer
if IsSkimTree:
  if not IsTAMSA:
    print "Skimming only possible in SNU"
    exit()

## Machine-dependent variables

if IsKNU:
  args.Queue = "cms"

## Make Sample List

InputSample_Data = ["DoubleMuon", "DoubleEG", "SingleMuon", "SingleElectron", "SinglePhoton", "MuonEG", "EGamma"]
AvailableDataPeriods = []
if args.Era == "2016preVFP":
  AvailableDataPeriods = ["B_ver2","C","D","E","F"]
elif args.Era == "2016postVFP":
  AvailableDataPeriods = ["F","G","H"]
elif args.Era == "2017":
  AvailableDataPeriods = ["B","C","D","E","F"]
elif args.Era == "2018":
  AvailableDataPeriods = ["A", "B","C","D"]
else:
  print "[SKFlat.py] Wrong Era : "+args.Era
  exit(1)

InputSamples = []
StringForHash = ""

## When using txt file for input (i.e., -l option)

if args.InputSampleList is not "":
  lines = open(args.InputSampleList)
  for line in lines:
    if "#" in line:
      continue
    line = line.strip('\n')
    InputSamples.append(line)
    StringForHash += line
else:
  if args.InputSample in InputSample_Data:
    if args.DataPeriod=="ALL":
      for period in AvailableDataPeriods:
        InputSamples.append(args.InputSample+":"+period)
        StringForHash += args.InputSample+":"+period
    elif args.DataPeriod in AvailableDataPeriods:
      InputSamples.append(args.InputSample+":"+args.DataPeriod)
      StringForHash += args.InputSample+":"+args.DataPeriod
  else:
    InputSamples.append(args.InputSample)
    StringForHash += args.InputSample
FileRangesForEachSample = []

## add flags to hash
for flag in Userflags:
  StringForHash += flag

## Get Random Number for webdir

random.seed(hash(StringForHash+timestamp+args.Era))
RandomNumber = int(random.random()*1000000)
str_RandomNumber = str(RandomNumber)
webdirname = timestamp+"_"+str_RandomNumber
webdirpathbase = SKFlatRunlogDir+'/www/SKFlatAnalyzerJobLogs/'+webdirname
while os.path.isdir(webdirpathbase):
  webdirpathbase += '_'

## skim string

SkimString = ""
if args.Skim!="":
  SkimString = args.Skim+"_"

## Define MasterJobDir

MasterJobDir = SKFlatRunlogDir+'/'+timestamp+'__'+str_RandomNumber+"__"+args.Analyzer+'__'+'Era'+args.Era
if args.Skim!="":
  MasterJobDir += "__"+args.Skim
for flag in Userflags:
  MasterJobDir += '__'+flag
MasterJobDir += '__'+HOSTNAME+'/'

## Copy libray

os.system('mkdir -p '+MasterJobDir+'/lib/')
os.system('cp '+SKFlat_LIB_PATH+'/* '+MasterJobDir+'/lib')

## Loop over samples

# true or false for each sample
SampleFinishedForEachSample = []
PostJobFinishedForEachSample = []
BaseDirForEachSample = []
XsecForEachSample = []
for InputSample in InputSamples:

  NJobs = args.NJobs

  SampleFinishedForEachSample.append(False)
  PostJobFinishedForEachSample.append(False)

  ## Global Varialbes

  IsDATA = False
  DataPeriod = ""
  if ":" in InputSample:
    IsDATA = True
    tmp = InputSample
    InputSample = tmp.split(":")[0]
    DataPeriod = tmp.split(":")[1]

  ## Prepare output

  base_rundir = MasterJobDir+InputSample
  if IsDATA:
    base_rundir = base_rundir+'_period'+DataPeriod
  base_rundir = base_rundir+"/"

  os.system('mkdir -p '+base_rundir)
  os.system('mkdir -p '+base_rundir+'/output/')

  ## Create webdir
  ## cf) base_rundir = $SKFlatRunlogDir/2019_02_26_222038__GetEffLumi__Era2016__KISTI/WW_pythia/

  this_webdir = webdirpathbase+'/'+base_rundir.replace(SKFlatRunlogDir,'').replace(HOSTNAME+'/',HOSTNAME+'__')
  os.system('mkdir -p '+this_webdir)

  ## If KNU, copy grid cert

  if IsKNU:
    os.system('cp /tmp/x509up_u'+UID+' '+base_rundir)

  ## Get Sample Path

  lines_files = []

  tmpfilepath = SAMPLE_DATA_DIR+'/For'+SampleHOSTNAME+'/'+SkimString+InputSample+'.txt'
  if IsDATA:
    tmpfilepath = SAMPLE_DATA_DIR+'/For'+SampleHOSTNAME+'/'+SkimString+InputSample+'_'+DataPeriod+'.txt'
  lines_files = os.popen("sed 's/#.*//' "+tmpfilepath+"|grep '.root'").readlines()
  os.system('cp '+tmpfilepath+' '+base_rundir+'/input_filelist.txt')

  NTotalFiles = len(lines_files)

  if NJobs>NTotalFiles or NJobs==0:
    NJobs = NTotalFiles

  SubmitOutput = open(base_rundir+'/SubmitOutput.log','w')

  SubmitOutput.write("<SKFlat> NTotalFiles = "+str(NTotalFiles)+'\n')
  SubmitOutput.write("<SKFlat> NJobs = "+str(NJobs)+'\n')
  nfilepjob = int(NTotalFiles/NJobs)
  SubmitOutput.write("<SKFlat> --> # of files per job = "+str(nfilepjob)+'\n')
  nfilepjob_remainder = NTotalFiles-(NJobs)*(nfilepjob)
  if nfilepjob_remainder>=(NJobs):
    SubmitOutput.write('nfilepjob_remainder = '+str(nfilepjob_remainder)+'\n')
    SubmitOutput.write('while, NJobs = '+str(NJobs)+'\n')
    SubmitOutput.write('--> exit'+'\n')
    sys.exit()

  FileRanges = []
  temp_end_largerjob = 0
  nfile_checksum = 0

  ## First nfilepjob_remainder jobs will have (nfilepjob+1) files per job

  for it_job in range(0,nfilepjob_remainder):
    FileRanges.append(range(it_job*(nfilepjob+1),(it_job+1)*(nfilepjob+1)))
    temp_end_largerjob = (it_job+1)*(nfilepjob+1)
    nfile_checksum += len(range(it_job*(nfilepjob+1),(it_job+1)*(nfilepjob+1)))

  ## Remaining NJobs-nfilepjob_remainder jobs will have (nfilepjob) files per job

  for it_job in range(0,NJobs-nfilepjob_remainder):
    FileRanges.append(range(temp_end_largerjob+(it_job*nfilepjob),temp_end_largerjob+((it_job+1)*nfilepjob) ))
    nfile_checksum += len(range(temp_end_largerjob+(it_job*nfilepjob),temp_end_largerjob+((it_job+1)*nfilepjob) ))
  SubmitOutput.write('nfile_checksum = '+str(nfile_checksum)+'\n')
  SubmitOutput.write('NTotalFiles = '+str(NTotalFiles)+'\n')
  SubmitOutput.close()
  FileRangesForEachSample.append(FileRanges)

  ## Get xsec and SumW

  this_dasname = ""
  this_xsec = -1
  this_sumw = -1
  if not IsDATA:
    if not os.path.exists(SAMPLE_DATA_DIR+'/CommonSampleInfo/'+InputSample+'.txt'):
      print("No "+SAMPLE_DATA_DIR+'/CommonSampleInfo/'+InputSample+'.txt')
      exit(2)
    lines_SamplePath = open(SAMPLE_DATA_DIR+'/CommonSampleInfo/'+InputSample+'.txt').readlines()
    for line in lines_SamplePath:
      if line[0]=="#":
        continue
      words = line.split()
      if InputSample==words[0]:
        this_dasname = words[1]
        this_xsec = words[2]
        this_sumw = words[4]
        break

  XsecForEachSample.append(this_xsec)

  ## Write run script

  if IsKISTI or IsTAMSA:

    commandsfilename = args.Analyzer+'_'+args.Era+'_'+InputSample
    if IsDATA:
      commandsfilename += '_'+DataPeriod
    for flag in Userflags:
      commandsfilename += '__'+flag
    run_commands = open(base_rundir+'/'+commandsfilename+'.sh','w')
    print>>run_commands,'''#!/bin/bash
SECTION=`printf $1`
WORKDIR=`pwd`

Trial=0

#### make sure use C locale
export LC_ALL=C

#### Don't make root history
export ROOT_HIST=0

#### use cvmfs for root ####
export CMS_PATH=/cvmfs/cms.cern.ch
source $CMS_PATH/cmsset_default.sh
export SCRAM_ARCH={2}
export cmsswrel={3}
cd /cvmfs/cms.cern.ch/$SCRAM_ARCH/cms/$cmsswrel/src
echo "@@@@ SCRAM_ARCH = "$SCRAM_ARCH
echo "@@@@ cmsswrel = "$cmsswrel
echo "@@@@ scram..."
eval `scramv1 runtime -sh`
cd -
source /cvmfs/cms.cern.ch/$SCRAM_ARCH/cms/$cmsswrel/external/$SCRAM_ARCH/bin/thisroot.sh

### modifying LD_LIBRARY_PATH to use libraries in base_rundir
export LD_LIBRARY_PATH=$(echo $LD_LIBRARY_PATH|sed 's@'$SKFlat_WD'/lib@{0}/lib@')

while [ "$Trial" -lt 3 ]; do
  echo "#### running ####"
  echo "root -l -b -q {1}/run_${{SECTION}}.C"
  root -l -b -q {1}/run_${{SECTION}}.C 2> err.log 
  EXITCODE=$?
  if [ "$EXITCODE" -eq 5 ]; then
    echo "IO error occured.. running again in 300 seconds.."
    Trial=$((Trial+=1))
    sleep 300
  else
    break
  fi
done

if [ "$EXITCODE" -ne 0 ]; then
  echo "ERROR errno=$EXITCODE" >> err.log
fi

cat err.log >&2
exit $EXITCODE
'''.format(MasterJobDir, base_rundir, SCRAM_ARCH, cmsswrel)
    run_commands.close()

    submit_command = open(base_rundir+'/submit.jds','w')
    if IsUI10:
      print>>submit_command,'''executable = {1}.sh
universe   = vanilla
arguments  = $(Process)
requirements = OpSysMajorVer == 6
log = condor.log
getenv     = True
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
output = job_$(Process).log
error = job_$(Process).err
transfer_output_remaps = "hists.root = output/hists_$(Process).root"
queue {0}
'''.format(str(NJobs), commandsfilename)
      submit_command.close()
    elif IsUI20:
      print>>submit_command,'''executable = {1}.sh
universe   = vanilla
requirements = ( HasSingularity == true )
arguments  = $(Process)
log = condor.log
getenv     = True
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
output = job_$(Process).log
error = job_$(Process).err
accounting_group=group_cms
+SingularityImage = "/cvmfs/singularity.opensciencegrid.org/opensciencegrid/osgvo-el6:latest"
+SingularityBind = "/cvmfs, /cms, /share"
transfer_output_remaps = "hists.root = output/hists_$(Process).root"
queue {0}
'''.format(str(NJobs), commandsfilename)
      submit_command.close()
    elif IsTAMSA:
      concurrency_limits=''
      if args.NMax:
        concurrency_limits='concurrency_limits = n'+str(args.NMax)+'.'+os.getenv("USER")
      request_memory=''
      if args.Memory:
        request_memory='request_memory = '+str(args.Memory)
      print>>submit_command,'''executable = {1}.sh
universe   = vanilla
arguments  = $(Process)
log = condor.log
getenv     = True
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
output = job_$(Process).log
error = job_$(Process).err
transfer_output_remaps = "hists.root = output/hists_$(Process).root"
{2}
{3}
queue {0}
'''.format(str(NJobs), commandsfilename,concurrency_limits,request_memory)
      submit_command.close()

  CheckTotalNFile=0
  for it_job in range(0,len(FileRanges)):

    #print "["+str(it_job)+"th]",
    #print FileRanges[it_job],
    #print " --> "+str(len(FileRanges[it_job]))

    CheckTotalNFile = CheckTotalNFile+len(FileRanges[it_job])

    thisjob_dir = base_rundir+'/job_'+str(it_job)+'/'

    runfunctionname = "run"
    libdir = (MasterJobDir+'/lib').replace('///','/').replace('//','/')+'/'
    runCfileFullPath = ""
    if IsKISTI or IsTAMSA:
      runfunctionname = "run_"+str(it_job)
      runCfileFullPath = base_rundir+'/run_'+str(it_job)+'.C'
    else:
      os.system('mkdir -p '+thisjob_dir)
      runCfileFullPath = thisjob_dir+'run.C'

    IncludeLine = 'R__LOAD_LIBRARY(/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/lhapdf/6.2.1-gnimlf3/lib/libLHAPDF.so)\n'

    out = open(runCfileFullPath, 'w')
    print>>out,'''{3}

void {2}(){{

  {0} m;

  m.SetTreeName("recoTree/SKFlat");
'''.format(args.Analyzer, libdir, runfunctionname, IncludeLine)

    out.write('  m.LogEvery = '+str(LogEvery)+';\n')

    if IsDATA:
      out.write('  m.IsDATA = true;\n')
      out.write('  m.DataStream = "'+InputSample+'";\n')
    else:
      out.write('  m.MCSample = "'+InputSample+'";\n');
      out.write('  m.IsDATA = false;\n')
      out.write('  m.xsec = '+str(this_xsec)+';\n')
      out.write('  m.sumW = '+str(this_sumw)+';\n')

      if args.FastSim:
        out.write('  m.IsFastSim = true;\n')
      else:
        out.write('  m.IsFastSim = false;\n')

    out.write('  m.SetEra("'+str(args.Era)+'");\n')

    if len(Userflags)>0:
      out.write('  m.Userflags = {\n')
      for flag in Userflags:
        out.write('    "'+flag+'",\n')
      out.write('  };\n')

    for it_file in FileRanges[it_job]:
      thisfilename = lines_files[it_file].strip('\n')
      out.write('  if(!m.AddFile("'+thisfilename+'")) exit(EIO);\n')

    if IsSkimTree:
      tmp_filename = lines_files[ FileRanges[it_job][0] ].strip('\n')
      ## /data7/DATA/SKFlat/v949cand2_2/2017/DATA/SingleMuon/periodB/181107_231447/0000/SKFlatNtuple_2017_DATA_100.root
      ## /data7/DATA/SKFlat/v949cand2_2/2017/MC/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/181108_152345/0000/SKFlatNtuple_2017_MC_100.root
      skimoutdir = '/gv0/DATA/SKFlat/'+SKFlatV+'/'+args.Era+'/'
      skimoutfilename = ""
      if IsDATA:
        skimoutdir += "DATA_"+args.Analyzer+"/"+InputSample+"/period"+DataPeriod+"/"
        skimoutfilename = "SKFlatNtuple_"+args.Era+"_DATA_"+str(it_job)+".root"
      else:
        skimoutdir += "MC_"+args.Analyzer+"/"+this_dasname+"/"
        skimoutfilename = "SKFlatNtuple_"+args.Era+"_MC_"+str(it_job)+".root"
      skimoutdir += timestamp+"/"

      os.system('mkdir -p '+skimoutdir)
      out.write('  m.SetOutfilePath("'+skimoutdir+skimoutfilename+'");\n')

    else:
      if IsKISTI or IsTAMSA:
        out.write('  m.SetOutfilePath("hists.root");\n')
      else:
        out.write('  m.SetOutfilePath("'+thisjob_dir+'/hists.root");\n')

    if args.Reduction>1:
      out.write('  m.MaxEvent=m.fChain->GetEntries()/'+str(args.Reduction)+';\n')

    print>>out,'''  m.Init();
  m.initializeAnalyzer();
  m.initializeAnalyzerTools();
  m.SwitchToTempDir();
  m.Loop();

  m.WriteHist();

}'''

    out.close()

    if IsKNU:
      run_commands = open(thisjob_dir+'commands.sh','w')
      print>>run_commands,'''cd {0}
cp ../x509up_u{1} /tmp/
echo "[SKFlat.py] Okay, let's run the analysis"
root -l -b -q run.C 1>stdout.log 2>stderr.log
'''.format(thisjob_dir,UID)
      run_commands.close()

      jobname = 'job_'+str(it_job)+'_'+args.Analyzer
      cmd = 'qsub -V -q '+args.Queue+' -N '+jobname+' commands.sh'

      if not args.no_exec:
        cwd = os.getcwd()
        os.chdir(thisjob_dir)
        os.system(cmd+' > submitlog.log')
        os.chdir(cwd)
      sublog = open(thisjob_dir+'/submitlog.log','a')
      sublog.write('\nSubmission command was : '+cmd+'\n')
      sublog.close()

  if IsKISTI or IsTAMSA:

    cwd = os.getcwd()
    os.chdir(base_rundir)
    if not args.no_exec:
      condorOptions = ''
      if args.BatchName!="":
        condorOptions = ' -batch-name '+args.BatchName
      os.system('condor_submit submit.jds '+condorOptions)
    os.chdir(cwd)

  else:

    if args.no_exec:
      continue

    ## Write Kill Command

    KillCommand = open(base_rundir+'/Script_JobKill.sh','w')
    for it_job in range(0,len(FileRanges)):
      thisjob_dir = base_rundir+'/job_'+str(it_job)+'/'
      jobid = GetJobID(thisjob_dir, args.Analyzer, it_job, HOSTNAME)
      KillCommand.write('qdel '+jobid+' ## job_'+str(it_job)+' ##\n')
    KillCommand.close()

if args.no_exec:
  print '- RunDir = '+base_rundir
  exit()

## Set Output directory
## if args.Outputdir is not set, go to default setting

FinalOutputPath = args.Outputdir
if args.Outputdir=="":
  FinalOutputPath = SKFlatOutputDir+'/'+SKFlatV+'/'+args.Analyzer+'/'+args.Era+'/'
  for flag in Userflags:
    FinalOutputPath += flag+"__"
  if IsDATA:
    FinalOutputPath += '/DATA/'
if IsSkimTree:
  FinalOutputPath = '/gv0/DATA/SKFlat/'+SKFlatV+'/'+args.Era+'/'

os.system('mkdir -p '+FinalOutputPath)

print '##################################################'
print 'Submission Finished'
print '- JobID = '+str_RandomNumber
print '- Analyzer = '+args.Analyzer
print '- Skim = '+args.Skim
print '- InputSamples =',
print InputSamples
print '- NJobs = '+str(NJobs)
print '- Era = '+args.Era
print '- UserFlags =',
print Userflags
if IsKNU:
  print '- Queue = '+args.Queue
print '- RunDir = '+base_rundir
print '- output will be send to : '+FinalOutputPath
print '##################################################'

##########################
## Submittion all done. ##
## Now monitor job      ##
##########################

## Loop over samples again

AllSampleFinished = False
GotError = False
ErrorLog = ""

try:
  while not AllSampleFinished:

    if GotError:
      break

    AllSampleFinished = True

    for it_sample in range(0,len(InputSamples)):

      InputSample = InputSamples[it_sample]
      SampleFinished = SampleFinishedForEachSample[it_sample]
      PostJobFinished = PostJobFinishedForEachSample[it_sample]

      if PostJobFinished:
        continue
      else:
        AllSampleFinished = False

      ## Global Varialbes

      IsDATA = False
      DataPeriod = ""
      if ":" in InputSample:
        IsDATA = True
        tmp = InputSample
        InputSample = tmp.split(":")[0]
        DataPeriod = tmp.split(":")[1]

      ## Prepare output
      ## This should be copied from above

      base_rundir = MasterJobDir+InputSample
      if IsDATA:
        base_rundir = base_rundir+'_period'+DataPeriod
      base_rundir = base_rundir+"/"

      ## base_rundir = $SKFlatRunlogDir/2019_02_26_222038__GetEffLumi__Era2016__KISTI/WW_pythia/

      this_webdir = webdirpathbase+'/'+base_rundir.replace(SKFlatRunlogDir,'').replace(HOSTNAME+'/',HOSTNAME+'__')

      if not SampleFinished:

        ## This sample was not finished in the previous monitoring
        ## Monitor again this time

        ThisSampleFinished = True

        ## Write Job status until it's done

        statuslog = open(base_rundir+'/JobStatus.log','w')
        statuslog.write('Job submitted at '+string_JobStartTime+'\n')
        statuslog.write('JobNumber\t| Status\n')

        ToStatuslog = []
        n_eventran = 0
        finished = []
        EventDone = 0
        EventTotal = 0

        TotalEventRunTime = 0
        MaxTimeLeft = 0
        MaxEventRunTime = 0

        FileRanges = FileRangesForEachSample[it_sample]

        for it_job in range(0,len(FileRanges)):

          thisjob_dir = base_rundir+'/'
          if IsKISTI or IsTAMSA:
            thisjob_dir = base_rundir

          this_status = ""
          this_status = CheckJobStatus(thisjob_dir, args.Analyzer, it_job, HOSTNAME)

          if "ERROR" in this_status:
            GotError = True
            statuslog.write("#### ERROR OCCURED ####\n")
            statuslog.write(this_status+'\n')
            ErrorLog = this_status
            break

          if "FINISHED" not in this_status:
            ThisSampleFinished = False

          outlog = ""
          if "FINISHED" in this_status:
            finished.append("Finished")

            EventInfo = this_status.split()[1].split(':')

            this_EventDone = int(EventInfo[2])
            this_EventTotal = int(EventInfo[2])

            EventDone += this_EventDone
            EventTotal += this_EventTotal

            #### start
            line_EventRunTime = this_status.split()[2]+' '+this_status.split()[3]
            this_jobstarttime = GetDatetimeFromMyFormat(line_EventRunTime)
            #### end
            line_EventEndTime = this_status.split()[4]+' '+this_status.split()[5]
            this_jobendtime   = GetDatetimeFromMyFormat(line_EventEndTime)

            this_diff = this_jobendtime-this_jobstarttime
            this_EventRunTime = 86400*this_diff.days+this_diff.seconds

            this_TimePerEvent = float(this_EventRunTime)/float(this_EventDone)
            this_TimeLeft = (this_EventTotal-this_EventDone)*this_TimePerEvent

            TotalEventRunTime += this_EventRunTime
            MaxTimeLeft = max(MaxTimeLeft,this_TimeLeft)
            MaxEventRunTime = max(MaxEventRunTime,this_EventRunTime)

          elif "RUNNING" in this_status:
            outlog = str(it_job)+'\t| '+this_status.split()[1]+' %'

            EventInfo = this_status.split()[2].split(':')

            this_EventDone = int(EventInfo[1])
            this_EventTotal = int(EventInfo[2])

            EventDone += this_EventDone
            EventTotal += this_EventTotal

            line_EventRunTime = this_status.split()[3]+' '+this_status.split()[4]
            this_jobstarttime = GetDatetimeFromMyFormat(line_EventRunTime)
            this_diff = datetime.datetime.now()-this_jobstarttime
            this_EventRunTime = 86400*this_diff.days+this_diff.seconds

            if this_EventDone==0:
              this_EventDone = 1

            this_TimePerEvent = float(this_EventRunTime)/float(this_EventDone)
            this_TimeLeft = (this_EventTotal-this_EventDone)*this_TimePerEvent

            TotalEventRunTime += this_EventRunTime
            MaxTimeLeft = max(MaxTimeLeft,this_TimeLeft)
            MaxEventRunTime = max(MaxEventRunTime,this_EventRunTime)

            round_this_TimeLeft = round(this_TimeLeft,1)
            round_this_EventRunTime = round(this_EventRunTime,1)

            outlog += ' ('+str(round_this_EventRunTime)+' s ran, and '+str(round_this_TimeLeft)+' s left)'
            ToStatuslog.append(outlog)
            n_eventran += 1

          else:
            outlog = str(it_job)+'\t| '+this_status
            ToStatuslog.append(outlog)

          ##---- END it_job loop

        if GotError:

          ## When error occured, change both Finished/PostJob Flag to True

          SampleFinishedForEachSample[it_sample] = True
          PostJobFinishedForEachSample[it_sample] = True
          break

        for l in ToStatuslog:
          statuslog.write(l+'\n')
        statuslog.write('\n==============================================================\n')
        statuslog.write('HOSTNAME = '+HOSTNAME+'\n')
        statuslog.write('queue = '+args.Queue+'\n')
        statuslog.write(str(len(FileRanges))+' jobs submitted\n')
        statuslog.write(str(n_eventran)+' jobs are running\n')
        statuslog.write(str(len(finished))+' jobs are finished\n')

        ThisTime = datetime.datetime.now()
        string_ThisTime =  ThisTime.strftime('%Y-%m-%d %H:%M:%S')

        statuslog.write('XSEC = '+str(XsecForEachSample[it_sample])+'\n')
        statuslog.write('EventDone = '+str(EventDone)+'\n')
        statuslog.write('EventTotal = '+str(EventTotal)+'\n')
        statuslog.write('EventLeft = '+str(EventTotal-EventDone)+'\n')
        statuslog.write('TotalEventRunTime = '+str(TotalEventRunTime)+'\n')
        statuslog.write('MaxTimeLeft = '+str(MaxTimeLeft)+'\n')
        statuslog.write('MaxEventRunTime = '+str(MaxEventRunTime)+'\n')

        t_per_event = 1
        if EventDone is not 0:
          t_per_event = float(TotalEventRunTime)/float(EventDone)
        statuslog.write('t_per_event = '+str(t_per_event)+'\n')

        EstTime = ThisTime+datetime.timedelta(0, MaxTimeLeft)

        statuslog.write('Estimated Finishing Time : '+EstTime.strftime('%Y-%m-%d %H:%M:%S')+'\n')
        statuslog.write('Last checked at '+string_ThisTime+'\n')
        statuslog.close()

        ## copy statuslog to webdir

        os.system('cp '+base_rundir+'/JobStatus.log '+this_webdir)

        ## This time, it is found to be finished
        ## Change the flag

        if ThisSampleFinished:
          SampleFinishedForEachSample[it_sample] = True

        ##---- END if finished

      else:

        ## Job was finished in the previous monitoring
        ## Check if PostJob is also finished

        if not PostJobFinished:

          ## PostJob was not done in the previous monitoring
          ## Copy output, and change the PostJob flag


          ## if Skim, no need to hadd. move on!

          if IsSkimTree:
            PostJobFinishedForEachSample[it_sample] = True
            continue

          outputname = args.Analyzer+'_'+SkimString+InputSample
          if IsDATA:
            outputname += '_'+DataPeriod

          if not GotError:
            cwd = os.getcwd()
            os.chdir(base_rundir)

            #### if number of job is 1, we can just move the file, not hadd
            nFiles = len( FileRangesForEachSample[it_sample] )
            if nFiles==1:
              if IsKISTI or IsTAMSA:
                os.system('echo "nFiles = 1, so skipping hadd and just move the file" >> JobStatus.log')
                os.system('ls -1 output/*.root >> JobStatus.log')
                os.system('mv output/hists_0.root '+outputname+'.root')
              else:
                os.system('echo "nFiles = 1, so skipping hadd and just move the file" >> JobStatus.log')
                os.system('ls -1 job_0/*.root >> JobStatus.log')
                os.system('mv job_0/hists.root '+outputname+'.root')

            else:
              if IsKISTI or IsTAMSA:
                os.system('hadd -f '+outputname+'.root output/*.root >> JobStatus.log')
                os.system('rm output/*.root')
                #os.system('condor_run -a request_cpus=10 "hadd -j 10 -f '+outputname+'.root output/*.root 2>&1 >> JobStatus.log"')
              else:
                os.system('hadd -f '+outputname+'.root job_*/*.root >> JobStatus.log')
                os.system('rm job_*/*.root')

            ## Final Outputpath

            os.system('mv '+outputname+'.root '+FinalOutputPath)
            os.chdir(cwd)

          PostJobFinishedForEachSample[it_sample] = True

    if SendLogToWeb:

      os.system('cp -r '+webdirpathbase+'/* '+SKFlatLogWebDir)

    time.sleep(20)

except KeyboardInterrupt:
  print('interrupted!')

## Send Email now

from SendEmail import *
JobFinishEmail = '''#### Job Info ####
HOST = {3}
JobID = {6}
Analyzer = {0}
Era = {7}
Skim = {5}
# of Jobs = {4}
InputSample = {1}
{8}
Output sent to : {2}
'''.format(args.Analyzer,InputSamples,FinalOutputPath,HOSTNAME,NJobs,args.Skim,str_RandomNumber,args.Era,GetXSECTable(InputSamples,XsecForEachSample))
JobFinishEmail += '''##################
Job started at {0}
Job finished at {1}
'''.format(string_JobStartTime,string_ThisTime)

if IsKNU:
  JobFinishEmail += 'Queue = '+args.Queue+'\n'

EmailTitle = '['+HOSTNAME+']'+' Summary of JobID '+str_RandomNumber
if GotError:
  JobFinishEmail = "#### ERROR OCCURED ####\n"+JobFinishEmail
  JobFinishEmail = ErrorLog+"\n------------------------------------------------\n"+JobFinishEmail
  EmailTitle = '[ERROR] Summary of JobID '+str_RandomNumber

if IsKNU:
  SendEmailbyGMail(USER,SKFlatLogEmail,EmailTitle,JobFinishEmail)
else:
  SendEmail(USER,SKFlatLogEmail,EmailTitle,JobFinishEmail)
