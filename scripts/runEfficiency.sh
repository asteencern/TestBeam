#!/bin/bash

run=$1
configuration=$2
ped=1200

path="/afs/cern.ch/user/a/asteen/cmssw/CMSSW_8_0_1/src/HGCal"
CFG="${path}/myTest_cfg.py"
PEDHIGHGAIN="${path}/CondObjects/data/pedHighGain"${ped}".txt"
PEDLOWGAIN="${path}/CondObjects/data/pedLowGain"${ped}".txt"
PYTHONSCRIPT="${path}/runNewEB.py"
BAD_SPILL_SCRIPT="${path}/scripts/badSpillHelper.py"
TOP=${PWD}

echo ${path}
echo ${CFG}
echo ${PEDHIGHGAIN}
echo ${PEDLOWGAIN}
echo ${PYTHONSCRIPT}
echo ${BAD_SPILL_SCRIPT}

cp ${CFG} .
cp ${PEDHIGHGAIN} .
cp ${PEDLOWGAIN} .
cp ${PYTHONSCRIPT} pythonscript.py
cp ${BAD_SPILL_SCRIPT} .

cd ${path}
eval `scramv1 runtime -sh` 
cd ${TOP}

python pythonscript.py --runType=HGCRun --process=efficiency --runNumber=${run} --nSpills=15 --configuration=${configuration} --pedestalRun=${ped} --pedestalPath=./
