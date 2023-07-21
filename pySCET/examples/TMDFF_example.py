import pySCET.Fragmentation_Function.OPEmatching as OPE
import numpy as np
import pySCET.Constants_and_functions.constants as C

nnff = OPE.TMDmatching(OPEtype='FF',initPDF=False)
nnff.setxQtable(np.array([0.1,0.5,1]),np.array([1,5,10]))
nnff.setpartontype('q')
nnff.setquarkflavourtable([-3,-2,-1,1,2,3])
nnff.setPDF(pdfname='NNFF10_PIsum_nlo',setnumber=0)
Cqq=lambda z:np.array([0,2*C.CF*(1-z),2*C.CF*(2*(1+z**2)*np.log(z))])
Cgq=lambda z: 2*C.CF*np.array([0,2*(1+(1-z)**2)/z*np.log(z)+z])
nnff.setMatchingCoeffecient(Cqq=Cqq,Cgq=Cgq)
nnff.ParallelRun()
nnff.output()

