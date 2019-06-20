#! /usr/bin/env python
import os
import glob
import math
from array import array
import sys
import time
import subprocess
import tarfile
import datetime
import commands

currentDir = os.getcwd();
CMSSWDir =  currentDir+"/../";

#inputFolder = "/store/user/lpcbacon/15/";
inputFolder = "/store/user/arapyan/Run2/";
TestRun = 0

doMC = True;
doData = False;
category = ["mu"];
#category = ["el","mu"];

lumi = 35900.0

changes = raw_input("\n\nWrite change summary: ")

print "==> ",changes

# Get date and time for output directory
## ADD "test" IN OUTPUT FOLDER IF YOU ARE TESTING SO THAT LATER YOU REMEMBER TO WHICH DIRECTORY YOU HAVE TO REMOVE FROM EOS
if TestRun:
	outputFolder = "/store/user/gchaudha/WZTree_"+datetime.datetime.now().strftime('%Y_%m_%d_%Hh%M')+"_TEST/";
	OutputLogPath = "OutPut_Logs/Logs_" + datetime.datetime.now().strftime('%Y_%m_%d_%Hh%M') + "_TEST";
else:
	outputFolder = "/store/user/gchaudha/WZTree_"+datetime.datetime.now().strftime('%Y_%m_%d_%Hh%M');
	OutputLogPath = "OutPut_Logs/Logs_" + datetime.datetime.now().strftime('%Y_%m_%d_%Hh%M');

print "Name of output dir: ",outputFolder
# create a directory on eos
os.system('xrdfs root://cmseos.fnal.gov/ mkdir ' + outputFolder)
# create directory in pwd for log files
os.system('mkdir -p ' + OutputLogPath + "/Logs")


# Function to create a tar file
def make_tarfile(output_filename, source_dir):
    with tarfile.open(output_filename, "w:gz") as tar:
            tar.add(source_dir, arcname=os.path.basename(source_dir))

# Get CMSSW directory path and name
cmsswDirPath = os.environ['CMSSW_BASE']
CMSSWRel = cmsswDirPath.split("/")[-1]
#cmsswDirPath = commands.getstatusoutput('echo ${CMSSW_BASE}')
#CMSSWRel = os.path.basename(cmsswDirPath[1])

print "CMSSW release used : ",CMSSWRel

os.system('echo "CMSSW Version used: '+CMSSWRel+'" > summary.dat')
os.system('echo "Current directory path: '+cmsswDirPath+'" >> summary.dat')
os.system('git diff > main.patch')
os.system("sed -i '1s/^/Changes Summary : "+changes+"\\n/' summary.dat")
os.system('echo -e "\n\n============\n== Latest commit summary \n\n" >> summary.dat ')
os.system("git log -1 --pretty=tformat:' Commit: %h %n Date: %ad %n Relative time: %ar %n Commit Message: %s' >> summary.dat ")
os.system('echo -e "\n\n============\n" >> summary.dat ')
os.system('git log -1 --format="%H" >> summary.dat ')
os.system('xrdcp -f summary.dat root://cmseos.fnal.gov/'+outputFolder+'/summary.dat')
os.system('xrdcp -f main.patch root://cmseos.fnal.gov/'+outputFolder+'/main.patch')

samples = [
 ( 0.0947916,	"WZTo3LNu_2Jets_MLL-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",       2280039, 0, "/Samples/190527_225923/0000/"),
 ( 0.1341072,	"WZTo3LNu_3Jets_MLL-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",      2080787, 0, "/Samples/190527_230429/0000/"),
 ( 0.4155876,	"WZTo3LNu_1Jets_MLL-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",      463529, 0, "/Samples/190527_225544/0000/"),
	#	Doubly Charged Higgs sample

	#	Singly Charged Higgs sample: Z -> LL

	#	Singly Charged Higgs sample: W -> Lv

    	#	EWK SM Signal

	#	aQGC Signal

    	##	QCD SM WWJJ
    	#( 5.568,	"WplusTo2JWminusToLNuJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	3994663,   0),

    	#	DY jets
    	##( 4435.5258,	"DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8",	0),
	#(0.0,	"DYJetsToLL_M_50_HT_100to200_TuneCP5_13TeV",	1,	0),
	#(0.0,	"DYJetsToLL_M_50_HT_200to400_TuneCP5_13TeV",	1,	0),
	#(0.0,	"DYJetsToLL_M_50_HT_400to600_TuneCP5_13TeV",	1,	0),
	#(0.0,	"DYJetsToLL_M_50_HT_600to800_TuneCP5_13TeV",	1,	0),
	#(0.0,	"DYJetsToLL_M_50_HT_800to1200_TuneCP5_13TeV",	1,	0),
	#(0.0,	"DYJetsToLL_M_50_HT_1200to2500_TuneCP5_13TeV",	1,	0),
	#(0.0,	"DYJetsToLL_M_50_HT_2500toInf_TuneCP5_13TeV",	1,	0),

    	##	Wjets
	#(0.0,	"WJetsToLNu_HT_100To200_TuneCP5_13TeV",		1,	0),
	#(0.0,	"WJetsToLNu_HT_200To400_TuneCP5_13TeV",		1,	0),
	#(0.0,	"WJetsToLNu_HT_400To600_TuneCP5_13TeV",		1,	0),
	#(0.0,	"WJetsToLNu_HT_600To800_TuneCP5_13TeV",		1,	0),
	#(0.0,	"WJetsToLNu_HT_800To1200_TuneCP5_13TeV",	1,	0),
	#(0.0,	"WJetsToLNu_HT_1200To2500_TuneCP5_13TeV",	1,	0),
	#(0.0,	"WJetsToLNu_HT_2500ToInf_TuneCP5_13TeV",	1,	0),

	##	VV/VVV
	##( 49.997,	"WWTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8",	5057358,	953706),
	#(0.0,	"WW_TuneCP5_13TeV_pythia8",	1,	0),
	#(0.0,	"WZ_TuneCP5_13TeV_pythia8",	1,	0),
	#(0.0,	"ZZ_TuneCP5_13TeV_pythia8",	1,	0),

	##	TTbar and single top
	#(0.0,	"ST_s_channel_4f_leptonDecays_TuneCP5_13TeV_amcatnlo_pythia8",				1,	0),
	#(0.0,	"ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV_powheg_pythia8",			1,	0),
	#(0.0,	"ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV_powheg_pythia8",				1,	0),
	#(0.0,	"ST_t_channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV_powhegV2_madspin_pythia8",	1,	0),
	#(0.0,	"ST_t_channel_top_4f_inclusiveDecays_TuneCP5_13TeV_powhegV2_madspin_pythia8",		1,	0),
	#(0.0,	"TTTo2L2Nu_TuneCP5_13TeV_powheg_pythia8",						1,	0),
	#(0.0,	"TTToHadronic_TuneCP5_13TeV_powheg_pythia8",						1,	0),
	#(0.0,	"TTToSemiLeptonic_TuneCP5_13TeV_powheg_pythia8",					1,	0),

	###	QCD
	##( 27990000.0,	"QCD_HT100to200_13TeV_1",	80614044,	1,	0),
	#(0.0,	"QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8",	1,	0),
	#(0.0,	"QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8",	1,	0),
	#(0.0,	"QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8",	1,	0),
	#(0.0,	"QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8",	1,	0),
	#(0.0,	"QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8",	1,	0),
    ]

nameDataMu = [
	"SingleMuonRun2017B_17Nov2017_v1",
	"SingleMuonRun2017C_17Nov2017_v1",
	"SingleMuonRun2017D_17Nov2017_v1",
	"SingleMuonRun2017E_17Nov2017_v1",
	"SingleMuonRun2017F_17Nov2017_v1"
   ];

nameDataEl = [
];


inputlist = "runstep2condor.sh, python/produceWWNtuples.py"

nameData = {"el": nameDataEl, "mu":nameDataMu};

command = "python python/produceWWNtuples.py -i "+inputFolder+" $*";

outScript = open("runstep2condor.sh","w");
outScript.write('#!/bin/bash');
outScript.write("\n"+'echo "Starting job on " `date`');
outScript.write("\n"+'echo "Running on: `uname -a`"');
outScript.write("\n"+'echo "System software: `cat /etc/redhat-release`"');
outScript.write("\n"+'source /cvmfs/cms.cern.ch/cmsset_default.sh');
outScript.write("\n"+'### copy the input root files if they are needed only if you require local reading');
outScript.write("\n"+'xrdcp -s root://cmseos.fnal.gov/'+outputFolder + '/'+CMSSWRel +'.tgz  .');
outScript.write("\n"+'tar -xf '+ CMSSWRel +'.tgz' );
outScript.write("\n"+'rm '+ CMSSWRel +'.tgz' );
outScript.write("\n"+'cd ' + CMSSWRel + '/src/WWAnalysis/WWAnalysisRun2' );
#outScript.write("\n"+'echo "====> List files : " ');
#outScript.write("\n"+'ls -alh');
outScript.write("\n"+'echo "====> Remove any file with name similar to WWTree*.root... " ');
outScript.write("\n"+'rm WWTree*.root');
outScript.write("\n"+'scramv1 b ProjectRename');
outScript.write("\n"+'eval `scram runtime -sh`');
#outScript.write("\n"+'echo "====> List files : " ');
#outScript.write("\n"+'ls -alh');
outScript.write("\n"+command);
#outScript.write("\n"+'echo "====> List files : " ');
#outScript.write("\n"+'ls -alh');
outScript.write("\n"+'echo "====> List root files : " ');
outScript.write("\n"+'ls WWTree*.root');
outScript.write("\n"+'echo "====> copying WWTree*.root file to stores area..." ');
outScript.write("\n"+'xrdcp -f WWTree*.root root://cmseos.fnal.gov/' + outputFolder);
#outScript.write("\n"+'xrdcp -f WWTree*.txt root://cmseos.fnal.gov/' + outputFolder);
outScript.write("\n"+'rm WWTree*.root');
outScript.write("\n"+'cd ${_CONDOR_SCRATCH_DIR}');
outScript.write("\n"+'rm -rf ' + CMSSWRel);
outScript.write("\n");
outScript.close();
os.system("chmod 777 runstep2condor.sh");

outJDL = open("runstep2condor.jdl","w");
outJDL.write("Executable = runstep2condor.sh\n");
outJDL.write("Universe = vanilla\n");
#outJDL.write("Requirements =FileSystemDomain==\"fnal.gov\" && Arch==\"X86_64\"");
outJDL.write("Notification = ERROR\n");
outJDL.write("Should_Transfer_Files = YES\n");
outJDL.write("WhenToTransferOutput = ON_EXIT\n");
#outJDL.write("include : list-infiles.sh |\n");
outJDL.write("Transfer_Input_Files = "+inputlist+"\n");
outJDL.write("x509userproxy = $ENV(X509_USER_PROXY)\n");


    #MC
if( doMC ):
    for i in range(len(samples)):
      #######################################################################
      ##
      ##		Text File Preparation	
      ##
      #######################################################################
      print("xrdfs root://cmseos.fnal.gov ls "+inputFolder+str(samples[i][1])+str(samples[i][4])+" | awk '{print \"root://cmseos.fnal.gov/\"$1}' > temp.txt")
      print "\n\n"
      os.system("xrdfs root://cmseos.fnal.gov ls "+inputFolder+str(samples[i][1])+str(samples[i][4])+" | awk '{print \"root://cmseos.fnal.gov/\"$1}' > temp.txt")
      tempFile = open('temp.txt','r')
      CountFiles = 0
      CountLines = 0
      for lines in tempFile:
        if CountLines % 10 == 0:
          newFile = open('listTemp_'+str(samples[i][1])+'_'+str(CountFiles)+'.txt','w')
          CountFiles = CountFiles + 1
        newFile.write(lines)
        CountLines = CountLines + 1
      ######################################################################
      #
      #		Main Part
      #
      ######################################################################
      for counts in xrange(CountFiles):
	#outJDL.write("Transfer_Input_Files = "+inputlist+" "+"listTemp_"+str(samples[i][1])+"_"+str(counts)+".txt  \n");
        outJDL.write("Output = "+OutputLogPath+"/"+str(samples[i][1])+"_"+str(counts)+".stdout\n");
        outJDL.write("Error  = "+OutputLogPath+"/"+str(samples[i][1])+"_"+str(counts)+".stdout\n");
        outJDL.write("Log  = "+OutputLogPath+"/Logs/"+str(samples[i][1])+"_"+str(counts)+".log\n");
        outJDL.write("Arguments = -n "+str(samples[i][1])+" -o WWTree_"+str(samples[i][1])+"_"+str(counts)+" -w "+str(samples[i][0])+" -no "+ str(samples[i][2]) + " -noNeg " + str(samples[i][3]) + " -lumi "+str(lumi)+" --ismc 1 -trig 1 -c lpc -f "+ "listTemp_"+str(samples[i][1])+"_"+str(counts)+".txt " +" \n");
        outJDL.write("Queue\n");
    
#data
if( doData ):
    for a in range(len(category)):
        for i in range(len(nameData[category[a]])):
      	    #######################################################################
      	    ##
      	    ##		Text File Preparation	
      	    ##
      	    #######################################################################
      	    print("xrdfs root://cmseos.fnal.gov ls "+inputFolder+str(samples[i][1])+" | awk '{print \"root://cmseos.fnal.gov/\"$1}' > temp.txt")
      	    print "\n\n"
      	    os.system("xrdfs root://cmseos.fnal.gov ls "+inputFolder+str(samples[i][1])+" | awk '{print \"root://cmseos.fnal.gov/\"$1}' > temp.txt")
      	    tempFile = open('temp.txt','r')
      	    CountFiles = 0
      	    CountLines = 0
      	    for lines in tempFile:
      	      if CountLines % 50 == 0:
      	        newFile = open('listTemp_'+str(samples[i][1])+'_'+str(CountFiles)+'.txt','w')
      	        CountFiles = CountFiles + 1
      	      newFile.write(lines)
      	      CountLines = CountLines + 1
      	    ######################################################################
      	    #
      	    #		Main Part
      	    #
      	    ######################################################################
      	    for counts in xrange(CountFiles):
               outJDL.write("Output = "+OutputLogPath+"/"+(nameData[category[a]])[i]+"_"+str(counts)+".stdout\n");
               outJDL.write("Error = "+OutputLogPath+"/"+(nameData[category[a]])[i]+"_"+str(counts)+".stdout\n");
               outJDL.write("Log = "+OutputLogPath+"/Logs/"+(nameData[category[a]])[i]+"_"+str(counts)+".log\n");
               outJDL.write("Arguments = -n "+(nameData[category[a]])[i]+" -o WWTree_"+(nameData[category[a]])[i]+"_"+category[a]+"_"+str(counts)+" -w 1. -no 1. --ismc 0 -trig 1 -c lpc -f "+ "listTemp_"+str(samples[i][1])+"_"+str(counts)+".txt " +" \n");
               print("Arguments = -n "+(nameData[category[a]])[i]+" -o WWTree_"+(nameData[category[a]])[i]+"_"+category[a]+"_"+str(counts)+" -w 1. -no 1. --ismc 0 -trig 1 -c lpc -f "+ "listTemp_"+str(samples[i][1])+"_"+str(counts)+".txt " +" \n");
               outJDL.write("Queue\n");

outJDL.close();

# create tarball of present working CMSSW base directory
os.system('rm -f CMSSW*.tgz')
make_tarfile(CMSSWRel+".tgz", cmsswDirPath)

# send the created tarball to eos
os.system('xrdcp -f ' + CMSSWRel+".tgz" + ' root://cmseos.fnal.gov/'+outputFolder+'/' + CMSSWRel+".tgz")

print "===> Delete the input text files"
#os.system("rm -f listTemp*.txt temp.txt")
print "===> Set Proxy Using:";
print "\tvoms-proxy-init --voms cms --valid 168:00";
print "\"condor_submit runstep2condor.jdl\" to submit";
