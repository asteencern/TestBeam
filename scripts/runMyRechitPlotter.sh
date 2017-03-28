#!/bin/bash

run=$1
ped=1200

path="/afs/cern.ch/user/a/asteen/cmssw/CMSSW_8_0_1/src/HGCal"
CFG="${path}/myTest_cfg.py"
PEDHIGHGAIN="${path}/CondObjects/data/pedHighGain"${ped}".txt"
PEDLOWGAIN="${path}/CondObjects/data/pedLowGain"${ped}".txt"
BAD_SPILL_SCRIPT="${path}/scripts/badSpillHelper.py"
PYTHONSCRIPT="${path}/runNewEB.py"
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

python pythonscript.py --runType=HGCRun --process=rechitplotter --runNumber=${run} --nSpills=15 --configuration=2 --pedestalRun=${ped} --pedestalPath=./
