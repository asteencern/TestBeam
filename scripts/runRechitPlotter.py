import os

#runs=[844, 844, 845, 846, 847, 848, 849, 850, 851, 852, 853, 854, 855, 856, 857, 858, 859, 860, 861, 862, 1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010, 1026, 1027, 1028, 1029, 1030, 1031, 1032, 1033, 1034, 1035, 1036, 1037, 1038, 1039, 1040, 1041, 1042, 1043, 1044, 1045, 1046, 1047]


#runs=[1312, 1313, 1314, 1315, 1316, 1317, 1318, 1319, 1320, 1321, 1322, 1323, 1324, 1325, 1326, 1327, 1328, 1329, 1330, 1331, 1332]

runs=[1411+i for i in range(0,242)] #1411->1652

for i in runs:
    print "bsub -q 1nh -o /afs/cern.ch/user/a/asteen/bsubOutput /afs/cern.ch/user/a/asteen/cmssw/CMSSW_8_0_1/src/HGCal/scripts/runRechitPlotter.sh "+str(i)
    os.system("bsub -q 1nh -o /afs/cern.ch/user/a/asteen/bsubOutput /afs/cern.ch/user/a/asteen/cmssw/CMSSW_8_0_1/src/HGCal/scripts/runRechitPlotter.sh "+str(i))
