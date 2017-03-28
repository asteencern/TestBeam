import os,csv

class badSpillFinder:
    filename=''

    def __init__(self,filename=''):
        self.filename=filename

    def Find(self):
        badspills=[]
        if self.filename=='':
            print "uncorrect file name"
            return
        with open(self.filename) as f:
            c = csv.reader(f, delimiter='\t', skipinitialspace=True)
            for line in c:
                findDanger=False
                if 'RUN' in line[0]:
                    for coll in line:
                        if 'DANGER=true' in coll:
                            findDanger=True
                            break
                    if findDanger==True:
                        badspill=int(line[1][6])+int(line[1][7])
                        if not badspill in badspills:
                            badspills.append(badspill)
        return badspills
    

class badSpillWriter:
    run=0
    badspills=[]

    def __init__(self,run,badspills):
        self.run=run
        self.badspills=badspills
        
    def Write(self):
        with open("bad_spill.txt",'w') as csvoutput:
            writer = csv.writer(csvoutput, delimiter='\t',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
            for i in self.badspills:
                writer.writerow([str(self.run),str(i)])
