import os
import util

#energies=[20,32,70,100,150,200,250]
#for i in energies:
#    u=util.util(i,"e-",True)
#    u.runMyRechitPlotter()



a=util.util(0,"mu-",True)
a.runMyRechitPlotter()

b=util.util(0,"mu-",False)
b.runMyRechitPlotter()

c=util.util(125,"pi-",True)
c.runMyRechitPlotter()
