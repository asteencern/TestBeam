import os,sys
from optparse import OptionParser

parser = OptionParser()

parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")

parser.add_option("-n", "--runNumber", dest="runNumber",type="int",
                  help="run number", default=1030)

parser.add_option("-t", "--runType", dest="runType",choices=["HGCRun","PED"],
                  help="run type",default="HGCRun")

parser.add_option("-p", "--process", dest="process",choices=["reco", "pedestal", "tracking", "rechitplotter","display","ntuple","efficiency"],
                  help="process to run",default="pedestal")

parser.add_option("-s", "--nSpills", dest="nSpills",type="int",
                  help="number of spill to run", default=3)

parser.add_option("-c", "--configuration", dest="configuration",type="int",
                  help=" hgcal setup configuration", default=1)

parser.add_option("-P", "--pedestalRun", dest="pedestalRun",type="int",
                  help="pedestal run number to use", default=1024)

parser.add_option("-C", "--pedestalPath", dest="pedestalPath",
                  help="path to find pedestal files", default="CondObjects/data/")

parser.add_option("-S", "--saveOnEos", dest="saveOnEos", type="int",
                  help="bool to set to save output on eps", default=1)

# parser.add_option("-e", "--nEvents", dest="nEvents",type="int",
#                   help="number of events to run", default=-1)

(options, args) = parser.parse_args()

print options

# example :
#
#  python runNewEB.py --runType=HGCRun --process=tracking --runNumber=1312 --nSpills=15 --configuration=2 --pedestalRun=1200 --pedestalPath=./
#
# 

cmd=""
eosSetup="/afs/cern.ch/project/eos/installation/cms/etc/setup.sh"
cmd+=eosSetup+"; "
eosPath="xrdcp -f root://eoscms.cern.ch//store/group/upgrade/HGCAL/TestBeam/CERN/Sept2016/"
cmd+=eosPath
eosFile="%s_Output_%06d.txt"%(options.runType,options.runNumber)
cmd+=eosFile+" "+eosFile

print eosPath, eosFile
os.system(cmd)

dataFolder='./'
outputFolder='./'
cmd=""
cmd="eval `scramv1 runtime -sh`;cmsRun myTest_cfg.py"
cmd+=" runNumber="+str(options.runNumber)+" nSpills="+str(options.nSpills)
cmd+=" dataFolder="+dataFolder
cmd+=" outputFolder="+outputFolder
cmd+=" runType="+options.runType
cmd+=" configuration="+str(options.configuration)
if options.process == "pedestal":
    pedestalsLowGain="./pedLowGain"+options.runNumber+".txt"
    pedestalsHighGain="./pedHighGain"+options.runNumber+".txt"
    cmd+=" pedestalsLowGain="+pedestalsLowGain
    cmd+=" pedestalsHighGain="+pedestalsHighGain
    cmd+=" chainSequence=1"
elif options.process == "rechitplotter":
    pedestalsLowGain=options.pedestalPath+"pedLowGain"+str(options.pedestalRun)+".txt"
    pedestalsHighGain=options.pedestalPath+"pedHighGain"+str(options.pedestalRun)+".txt"
    cmd+=" pedestalsLowGain="+pedestalsLowGain
    cmd+=" pedestalsHighGain="+pedestalsHighGain
    cmd+=" chainSequence=5"
elif options.process == "reco":
    pedestalsLowGain=options.pedestalPath+"pedLowGain"+str(options.pedestalRun)+".txt"
    pedestalsHighGain=options.pedestalPath+"pedHighGain"+str(options.pedestalRun)+".txt"
    cmd+=" pedestalsLowGain="+pedestalsLowGain
    cmd+=" pedestalsHighGain="+pedestalsHighGain
    cmd+=" chainSequence=4"
elif options.process == "tracking":
    pedestalsLowGain=options.pedestalPath+"pedLowGain"+str(options.pedestalRun)+".txt"
    pedestalsHighGain=options.pedestalPath+"pedHighGain"+str(options.pedestalRun)+".txt"
    cmd+=" pedestalsLowGain="+pedestalsLowGain
    cmd+=" pedestalsHighGain="+pedestalsHighGain
    cmd+=" chainSequence=8"
elif options.process == "display":
    pedestalsLowGain=options.pedestalPath+"pedLowGain"+str(options.pedestalRun)+".txt"
    pedestalsHighGain=options.pedestalPath+"pedHighGain"+str(options.pedestalRun)+".txt"
    cmd+=" pedestalsLowGain="+pedestalsLowGain
    cmd+=" pedestalsHighGain="+pedestalsHighGain
    cmd+=" chainSequence=3 maxEvents=10"
elif options.process == "ntuple":
    pedestalsLowGain=options.pedestalPath+"pedLowGain"+str(options.pedestalRun)+".txt"
    pedestalsHighGain=options.pedestalPath+"pedHighGain"+str(options.pedestalRun)+".txt"
    cmd+=" pedestalsLowGain="+pedestalsLowGain
    cmd+=" pedestalsHighGain="+pedestalsHighGain
    cmd+=" chainSequence=7"
elif options.process == "efficiency":
    pedestalsLowGain=options.pedestalPath+"pedLowGain"+str(options.pedestalRun)+".txt"
    pedestalsHighGain=options.pedestalPath+"pedHighGain"+str(options.pedestalRun)+".txt"
    cmd+=" pedestalsLowGain="+pedestalsLowGain
    cmd+=" pedestalsHighGain="+pedestalsHighGain
    cmd+=" chainSequence=9"

print cmd
os.system(cmd)
os.system("rm "+eosFile)

if options.process == "tracking" and options.saveOnEos==1:
    os.system("")
    eosOutputFile="%s_Output_%06d.root"%(options.runType,options.runNumber)
    os.system("source /afs/cern.ch/project/eos/installation/user/etc/setup.sh; xrdcp -f Output.root root://eosuser.cern.ch//eos/user/a/asteen/hgcal/data/sep2016/tracking-reco/"+eosOutputFile)
if options.process == "reco" and options.saveOnEos==1:
    os.system("")
    eosOutputFile="%s_Output_%06d.root"%(options.runType,options.runNumber)
    os.system("source /afs/cern.ch/project/eos/installation/user/etc/setup.sh; xrdcp -f Output.root root://eosuser.cern.ch//eos/user/a/asteen/hgcal/data/sep2016/shower-reco/"+eosOutputFile)
if options.process == "rechitplotter" and options.saveOnEos==1:
    os.system("")
    eosOutputFile="%s_Output_%06d.root"%(options.runType,options.runNumber)
    os.system("source /afs/cern.ch/project/eos/installation/user/etc/setup.sh; xrdcp -f Output.root root://eosuser.cern.ch//eos/user/a/asteen/hgcal/data/sep2016/rechitplotter/"+eosOutputFile)
if options.process == "efficiency" and options.saveOnEos==1:
    os.system("")
    eosOutputFile="%s_Output_%06d.root"%(options.runType,options.runNumber)
    os.system("source /afs/cern.ch/project/eos/installation/user/etc/setup.sh; xrdcp -f Output.root root://eosuser.cern.ch//eos/user/a/asteen/hgcal/data/sep2016/efficiency/"+eosOutputFile)
