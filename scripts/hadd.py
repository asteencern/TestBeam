import os
import util

os.system("source /afs/cern.ch/project/eos/installation/user/etc/setup.sh")

# e- data (ShowerAnalyzer)
#energies=[20,32,70,100,150,200,250]
#for i in energies:
#    u=util.util(i,"e-",True,1)
#    output="hgcal_e-_"+str(i)+"GeV_conf1.root"
#    eosPath="root://eosuser.cern.ch//eos/user/a/asteen/hgcal/data/sep2016/shower-reco"
#    outputPath="/afs/cern.ch/work/a/asteen/public/data/resultRoot/sep2016/shower-reco/"  
#    u.hadd(output,eosPath,outputPath)
#
#for i in energies:
#    u=util.util(i,"e-",True,0)
#    output="hgcal_e-_"+str(i)+"GeV_conf0.root"
#    eosPath="root://eosuser.cern.ch//eos/user/a/asteen/hgcal/data/sep2016/shower-reco"
#    outputPath="/afs/cern.ch/work/a/asteen/public/data/resultRoot/sep2016/shower-reco/"  
#    u.hadd(output,eosPath,outputPath)


# muon data (TrackAnalyzer)
#u=util.util(0,"mu-",True)
#output="hgcalMuonBeam.root"
#eosPath="root://eosuser.cern.ch//eos/user/a/asteen/hgcal/data/sep2016/tracking-reco"
#outputPath="/afs/cern.ch/work/a/asteen/public/data/resultRoot/sep2016/tracking/"  
#u.hadd(output,eosPath,outputPath)
#
#u=util.util(0,"mu-",False)
#output="hgcalMuonNoBeam.root"
#eosPath="root://eosuser.cern.ch//eos/user/a/asteen/hgcal/data/sep2016/tracking-reco"
#outputPath="/afs/cern.ch/work/a/asteen/public/data/resultRoot/sep2016/tracking/"  
#u.hadd(output,eosPath,outputPath)

# pion data (TrackAnalyzer)
#u=util.util(125,"pi-",True)
#output="hgcal125GeVPionBeam.root"
#eosPath="root://eosuser.cern.ch//eos/user/a/asteen/hgcal/data/sep2016/tracking-reco"
#outputPath="/afs/cern.ch/work/a/asteen/public/data/resultRoot/sep2016/tracking/"  
#u.hadd(output,eosPath,outputPath)
#

# # MyRechitPlotter
# # e- data
# energies=[20,32,70,100,150,200,250]
# for i in energies:
#     u=util.util(i,"e-",True)
#     output="hgcal_e-_"+str(i)+"GeV.root"
#     eosPath="root://eosuser.cern.ch//eos/user/a/asteen/hgcal/data/sep2016/rechitplotter"
#     outputPath="/afs/cern.ch/work/a/asteen/public/data/resultRoot/sep2016/rechitplotter/"  
#     u.hadd(output,eosPath,outputPath)


# muon data
u=util.util(0,"mu-",True)
output="hgcalMuonBeam.root"
eosPath="root://eosuser.cern.ch//eos/user/a/asteen/hgcal/data/sep2016/rechitplotter"
outputPath="/afs/cern.ch/work/a/asteen/public/data/resultRoot/sep2016/rechitplotter/"  
u.hadd(output,eosPath,outputPath)

u=util.util(0,"mu-",False)
output="hgcalMuonNoBeam.root"
eosPath="root://eosuser.cern.ch//eos/user/a/asteen/hgcal/data/sep2016/rechitplotter"
outputPath="/afs/cern.ch/work/a/asteen/public/data/resultRoot/sep2016/rechitplotter/"  
u.hadd(output,eosPath,outputPath)

# pion data
u=util.util(125,"pi-",True)
output="hgcal125GeVPionBeam.root"
eosPath="root://eosuser.cern.ch//eos/user/a/asteen/hgcal/data/sep2016/rechitplotter"
outputPath="/afs/cern.ch/work/a/asteen/public/data/resultRoot/sep2016/rechitplotter/"  
u.hadd(output,eosPath,outputPath)


# Cell efficiency analyzer
#eosPath="root://eosuser.cern.ch//eos/user/a/asteen/hgcal/data/sep2016/efficiency"
#outputPath="/afs/cern.ch/work/a/asteen/public/data/resultRoot/sep2016/efficiency/"  
#
## muon data
#u=util.util(0,"mu-",True)
#output="hgcalMuonBeam.root"
#u.hadd(output,eosPath,outputPath)
#
#u=util.util(0,"mu-",False)
#output="hgcalMuonNoBeam.root"
#u.hadd(output,eosPath,outputPath)

# pion data
#u=util.util(125,"pi-",True)
#output="hgcal125GeVPionBeam.root"
#u.hadd(output,eosPath,outputPath)

