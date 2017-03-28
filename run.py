import os,sys
import badSpillHelper
from optparse import OptionParser

parser = OptionParser()

parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")

parser.add_option("-n", "--runNumber", dest="runNumber",type="int",
                  help="run number", default=1030)

parser.add_option("-t", "--runType", dest="runType",choices=["HGCRun","PED"],
                  help="run type",default="HGCRun")

parser.add_option("-p", "--process", dest="process",choices=["reco", "pedestal", "tracking","display","ntuple"],
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

(options, args) = parser.parse_args()

print options

# example :
#
#  python run.py --runType=HGCRun --process=ntuple --runNumber=1312 --nSpills=15 --configuration=2 --pedestalRun=1200 --pedestalPath=./
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

finder=badSpillHelper.badSpillFinder(eosFile)
badspills=finder.Find()
writer=badSpillHelper.badSpillWriter(options.runNumber,badspills)
writer.Write()

print "Bad Spills = ", badspills

dataFolder='./'
outputFolder='./'
cmd=""
cmd="eval `scramv1 runtime -sh`;cmsRun cmssw_cfg.py"
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
elif options.process == "ntuple":
    pedestalsLowGain=options.pedestalPath+"pedLowGain"+str(options.pedestalRun)+".txt"
    pedestalsHighGain=options.pedestalPath+"pedHighGain"+str(options.pedestalRun)+".txt"
    cmd+=" pedestalsLowGain="+pedestalsLowGain
    cmd+=" pedestalsHighGain="+pedestalsHighGain
    cmd+=" chainSequence=2"
elif options.process == "display":
    pedestalsLowGain=options.pedestalPath+"pedLowGain"+str(options.pedestalRun)+".txt"
    pedestalsHighGain=options.pedestalPath+"pedHighGain"+str(options.pedestalRun)+".txt"
    cmd+=" pedestalsLowGain="+pedestalsLowGain
    cmd+=" pedestalsHighGain="+pedestalsHighGain
    cmd+=" chainSequence=3"
elif options.process == "tracking":
    pedestalsLowGain=options.pedestalPath+"pedLowGain"+str(options.pedestalRun)+".txt"
    pedestalsHighGain=options.pedestalPath+"pedHighGain"+str(options.pedestalRun)+".txt"
    cmd+=" pedestalsLowGain="+pedestalsLowGain
    cmd+=" pedestalsHighGain="+pedestalsHighGain
    cmd+=" chainSequence=6"
elif options.process == "reco":
    pedestalsLowGain=options.pedestalPath+"pedLowGain"+str(options.pedestalRun)+".txt"
    pedestalsHighGain=options.pedestalPath+"pedHighGain"+str(options.pedestalRun)+".txt"
    cmd+=" pedestalsLowGain="+pedestalsLowGain
    cmd+=" pedestalsHighGain="+pedestalsHighGain
    cmd+=" chainSequence=5"

print cmd
os.system(cmd)
