import os
import runList

class util:
    energy=0
    particle="e-"
    beam=True
    runlist=runList.runList()
    conf=1
    def __init__(self,energy=0,particle="e-",beam=True,conf=1):
        self.energy=energy
        self.particle=particle
        self.beam=beam
        self.conf=conf
        self.runlist=runList.runList(self.energy, self.particle, self.beam,self.conf)
        
    def printInfo(self):
        self.runlist.Print()

    def runShowerAnalysis(self):
        print "start shower analysis with : "
        self.printInfo()
        for i in self.runlist.runs:
            print "bsub -q 1nh -o /afs/cern.ch/user/a/asteen/bsubOutput /afs/cern.ch/user/a/asteen/cmssw/CMSSW_8_0_1/src/HGCal/scripts/runShower.sh "+str(i)+" "+str(self.conf+1)
            os.system("bsub -q 1nh -o /afs/cern.ch/user/a/asteen/bsubOutput /afs/cern.ch/user/a/asteen/cmssw/CMSSW_8_0_1/src/HGCal/scripts/runShower.sh "+str(i)+" "+str(self.conf+1))

    def runTrackAnalysis(self):
        print "start track analysis with : "
        self.printInfo()
        for i in self.runlist.runs:
            print "bsub -q 1nh -o /afs/cern.ch/user/a/asteen/bsubOutput /afs/cern.ch/user/a/asteen/cmssw/CMSSW_8_0_1/src/HGCal/scripts/runTracking.sh "+str(i)+" "+str(self.conf+1)
            os.system("bsub -q 1nh -o /afs/cern.ch/user/a/asteen/bsubOutput /afs/cern.ch/user/a/asteen/cmssw/CMSSW_8_0_1/src/HGCal/scripts/runTracking.sh "+str(i)+" "+str(self.conf+1))

    def runMyRechitPlotter(self):
        print "start shower analysis with : "
        self.printInfo()
        for i in self.runlist.runs:
            print "bsub -q 1nh -o /afs/cern.ch/user/a/asteen/bsubOutput /afs/cern.ch/user/a/asteen/cmssw/CMSSW_8_0_1/src/HGCal/scripts/runMyRechitPlotter.sh "+str(i)
            os.system("bsub -q 1nh -o /afs/cern.ch/user/a/asteen/bsubOutput /afs/cern.ch/user/a/asteen/cmssw/CMSSW_8_0_1/src/HGCal/scripts/runMyRechitPlotter.sh "+str(i))

    def hadd(self,outputName,eosPath,outputPath="./"):
        top=os.getcwd()
        os.system("source /afs/cern.ch/project/eos/installation/user/etc/setup.sh")
        os.chdir( outputPath )
        hadd="hadd -f "
        hadd+=outputName+" "
        rm="rm "
        for i in self.runlist.runs:
            eosFile="HGCRun_Output_%06d.root "%(i)
            cmd="xrdcp -f "+eosPath+"/"
            cmd+=eosFile+eosFile
            os.system( cmd )
            hadd+=eosFile
            rm+=eosFile

        print hadd
        os.system(hadd)
        os.system(rm)
        os.chdir( top )

    def runCellEfficiencyAnalysis(self):
        print "start efficiency analysis with : "
        self.printInfo()
        for i in self.runlist.runs:
            print "bsub -q 1nh -o /afs/cern.ch/user/a/asteen/bsubOutput /afs/cern.ch/user/a/asteen/cmssw/CMSSW_8_0_1/src/HGCal/scripts/runEfficiency.sh "+str(i)+" "+str(self.conf+1)
            os.system("bsub -q 1nh -o /afs/cern.ch/user/a/asteen/bsubOutput /afs/cern.ch/user/a/asteen/cmssw/CMSSW_8_0_1/src/HGCal/scripts/runEfficiency.sh "+str(i)+" "+str(self.conf+1))
