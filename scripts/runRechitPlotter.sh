#!/bin/bash

run=$1

path="/afs/cern.ch/user/a/asteen/cmssw/CMSSW_8_0_1/src/HGCal"
CFG="${path}/test_cfg.py"
PEDHIGHGAIN="${path}/CondObjects/data/pedHighGain1024.txt"
PEDLOWGAIN="${path}/CondObjects/data/pedLowGain1024.txt"
TOP=${PWD}

echo ${path}
echo ${CFG}
echo ${PEDHIGHGAIN}
echo ${PEDLOWGAIN}

cd ${path}
eval `scramv1 runtime -sh` 
cd ${TOP}
cmsRun ${CFG} runType=HGCRun pedestalsLowGain=${PEDLOWGAIN} pedestalsHighGain=${PEDHIGHGAIN} runNumber=$run nSpills=15 chainSequence=3 
