import os
import util

energies=[20,32,70,100,150,200,250]
for i in energies:
    u=util.util(i,"e-",True,1)
    u.runShowerAnalysis()

for i in energies:
    u=util.util(i,"e-",True,0)
    u.runShowerAnalysis()

