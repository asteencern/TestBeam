import os

class runList:
    energy=0
    particle="e-"
    beam=True
    conf=1
    runs=[]
    def __init__(self,energy=0,particle="e-",beam=True,conf=1):
        self.energy=energy
        self.particle=particle
        self.runs=[]
        self.beam=beam
        self.conf=conf
        # electron beam; conf 1 (up to 25 X0):
        if self.particle=="e-" and self.beam==True and self.conf==1:
            if self.energy==20:
                for i in range(1248,1271):
                    self.runs.append(i)

            elif self.energy==32:
                for i in range(1226,1248):
                    self.runs.append(i)
                #for i in range(1373,1386): #don't understand yet the problem
                #    self.runs.append(i)

            elif self.energy==70:
                for i in range(1153,1174):
                    self.runs.append(i)

            elif self.energy==100:
                for i in range(1202,1209):
                    self.runs.append(i)
                for i in range(1213,1224):
                    self.runs.append(i)
                #for i in range(1357,1373): #don't understand yet the problem
                #    self.runs.append(i)

            elif self.energy==150:
                for i in range(1122,1131):
                    self.runs.append(i)
                for i in range(1141,1152):
                    self.runs.append(i)

            elif self.energy==200:
                for i in range(1179,1200):
                    self.runs.append(i)
                #for i in range(1333,1345):
                #    self.runs.append(i)
                self.runs.remove(1197)
                self.runs.remove(1184)
            
            elif self.energy==250:
                for i in range(1291,1311):
                    self.runs.append(i)

        # electron beam; conf 0 (up to 15 X0):
        elif self.particle=="e-" and self.beam==True and self.conf==0:

            if self.energy==20:
                for i in range(916,935):
                    self.runs.append(i)

            elif self.energy==32:
                for i in range(1087,1109):
                    self.runs.append(i)
                self.runs.remove(1089)
                self.runs.remove(1090)

            elif self.energy==70:
                for i in range(866,890):
                    self.runs.append(i)
                self.runs.remove(870)
                self.runs.remove(871)
                self.runs.remove(884)
                self.runs.remove(885)
                self.runs.remove(886)

            elif self.energy==100:
                for i in range(1051,1074):
                    self.runs.append(i)
                self.runs.remove(1060)

            elif self.energy==150:
                for i in range(843,863):
                    self.runs.append(i)

            elif self.energy==200:
                for i in range(1026,1048):
                    self.runs.append(i)

            elif self.energy==250:
                for i in range(937,962):
                    self.runs.append(i)
                self.runs.remove(958)

        # muon runs
        elif self.particle=="mu-":
            if self.beam==True:
                for i in range(1312,1333):
                    self.runs.append(i)
            else:
                for i in range(1411,1652):
                    self.runs.append(i)
                self.runs.remove(1638)
                                                                
        # pion runs
        elif self.particle=="pi-" and self.beam==True:
            if self.energy==125:
                rmlist=[999,1012,1022]
                for i in range(998,1024):
                    if i not in rmlist:
                        self.runs.append(i)

    
    def addRun(self,run):
        self.runs.append(run)

    def removeRun(self,run):
        self.runs.remove(run)

    def Print(self):
        if self.energy!=0:
            print "hgcal",self.particle,"data at",self.energy," GeV"
        else:
            if self.beam==True:
                print "hgcal",self.particle,"data with beam"
            else:
                print "hgcal",self.particle,"data without beam"
        print "run list =", self.runs,"\n"
