
import lhapdf
import scipy.integrate as inte
import numpy as np
from p_tqdm import p_map
from itertools import chain



class OPE(object):
    def __init__(self):
        pass
    
    
    


def listtotitle(li):
    ti = ' '
    for i in li:
        ti += str(i) + ' '
    return ti





class TMDmatching(OPE):
    def __init__(self, OPEtype, initPDF = False):
        if OPEtype=='FF':
            self.OPEtype = OPEtype
            print('OPEtype set to be Fragmentation function')
        elif OPEtype=='PDF':
            self.OPEtype = OPEtype
            print('OPEtype set to be Parton Distribution function')
        else:
            print("OPEtype setting error, OPEtype should be FF or PDF. FF will be used as default, it will have no influence if collinear PDF(FF) and matching coeffecients are manually set.")
            self.OPEtype = 'FF'
        
        if initPDF:
            self.pdfname = "DSS07lo_pi"
            self.setnumber = 0
            self.pdf = lhapdf.mkPDF(self.pdfname, self.setnumber)
            
        self.x_value = [1.00000E-06,1.28121E-06,1.64152E-06,2.10317E-06,2.69463E-06,3.45242E-06,4.42329E-06,5.66715E-06,7.26076E-06,9.30241E-06,1.19180E-05,1.52689E-05,1.95617E-05,2.50609E-05,3.21053E-05,4.11287E-05,5.26863E-05,6.74889E-05,8.64459E-05,1.10720E-04,1.41800E-04,1.81585E-04,2.32503E-04,2.97652E-04,3.80981E-04,4.87518E-04,6.26039E-04,8.00452E-04,1.02297E-03,1.30657E-03,1.66759E-03,2.12729E-03,2.71054E-03,3.44865E-03,4.37927E-03,5.54908E-03,7.01192E-03,8.83064E-03,1.10763E-02,1.38266E-02,1.71641E-02,2.11717E-02,2.59364E-02,3.15062E-02,3.79623E-02,4.53425E-02,5.36750E-02,6.29705E-02,7.32221E-02,8.44039E-02,9.64793E-02,1.09332E-01,1.23067E-01,1.37507E-01,1.52639E-01,1.68416E-01,1.84794E-01,2.01731E-01,2.19016E-01,2.36948E-01,2.55242E-01,2.73927E-01,2.92954E-01,3.12340E-01,3.32036E-01,3.52019E-01,3.72282E-01,3.92772E-01,4.13533E-01,4.34326E-01,4.55495E-01,4.76836E-01,4.98342E-01,5.20006E-01,5.41818E-01,5.63773E-01,5.85861E-01,6.08077E-01,6.30459E-01,6.52800E-01,6.75387E-01,6.98063E-01,7.20830E-01,7.43683E-01,7.66623E-01,7.89636E-01,8.12791E-01,8.35940E-01,8.59175E-01,8.82485E-01,9.05866E-01,9.29311E-01,9.52817E-01,9.76387E-01,1.00000E+00]
        self.Q_value = [1.14018E+00,1.22539E+00,1.32482E+00,1.44156E+00,1.57954E+00,1.74381E+00,1.94091E+00,2.12132E+00,2.17945E+00,2.49622E+00,2.89383E+00,3.39925E+00,4.05063E+00,4.90286E+00,6.03623E+00,7.57063E+00,9.68869E+00,1.26749E+01,1.69835E+01,2.33578E+01,3.30501E+01,4.82334E+01,7.28041E+01,1.13998E+02,1.85779E+02,3.16228E+02]
        self.f_value = [-5,-4,-3,-2,-1,1,2,3,4,5,21]
        self.qf_value = [-5,-4,-3,-2,-1,1,2,3,4,5]
        self.gf_value = 21
        self.parton = 'qg'
        self.Cqq = lambda x: [0,0,0]
        self.Cgq = lambda x: [0,0,0]
        self.Cgg = lambda x: [0,0,0]
        self.Cqg = lambda x: [0,0,0]
        self.Cqpq = lambda x: [0,0,0]
        self.Cqbq = lambda x: [0,0,0]
        self.filename = 'TMD'+self.OPEtype+'.dat'
        
        
        
  
    
    def setPDF(self, pdfname, setnumber = 0):
        pdf = lhapdf.mkPDF(pdfname, setnumber)
        self.pdfname = pdfname
        self.setnumber = setnumber
        self.pdf = pdf
        self.filename = self.pdfname + '_TMD'+self.OPEtype+'_'+str(setnumber + 1000)+'.dat'
        
    
    
    def setMatchingCoeffecient(self, Cqq=lambda x: [0,0,0],Cgq=lambda x: [0,0,0],Cgg=lambda x: [0,0,0],Cqg=lambda x: [0,0,0],Cqpq=lambda x: [0,0,0],Cqbq=lambda x: [0,0,0]):
        self.Cqq = Cqq
        self.Cgq = Cgq
        self.Cgg = Cgg
        self.Cqg = Cqg
        self.Cqpq = Cqpq
        self.Cqbq = Cqbq
        # matching coeffecient should be an array being the (function) coeffecient of [delta(1-z), regular, plus(0,1-z),plus(1,1-z),...]
    

                        
    def matching(self,z,Q,flavour):
        if self.parton == 'q' and flavour == 21:
            return 0
        if self.parton == 'g' and flavour != 21:
            return 0
        if flavour == 21:
            gg = self.convolute(self.Cgg,lambda x: self.PDF(self.gf_value,x,Q),z)
            qg = sum([self.convolute(self.Cqg,lambda x: self.PDF(i,x,Q),z) for i in self.qf_value])
            return (gg+qg)*z
        else:
            gq = self.convolute(self.Cgq, lambda x: self.PDF(self.gf_value,x,Q),z)
            qq = self.convolute(self.Cqq, lambda x: self.PDF(flavour,x,Q),z)
            qbq = self.convolute(self.Cqbq, lambda x: self.PDF(-flavour,x,Q),z)
            qpq = 0
            for f in self.qf_value:
                if f!=flavour and f!=-flavour:
                    qpq +=self.convolute(self.Cqpq, lambda x: self.PDF(f,x,Q),z)
            return (gq+qq+qbq+qpq)*z
    
    ### *z is to match the convention of lhapdf
    
    def ParallelRun(self):
        f_value = self.f_value
        x_value = self.x_value
        Q_value = self.Q_value

        fcopy = f_value *len(x_value)*len(Q_value)
        b = [[Q_value[i]]*len(f_value) for i in range(len(Q_value))]
        qcopy=list(chain(*b)) * len(x_value)
        a = [[x_value[i]]*len(Q_value)*len(f_value) for i in range(len(x_value))]
        xcopy=list(chain(*a))
        
        parton = self.parton
        Cqq = self.Cqq
        Cqg=self.Cqg
        Cgq=self.Cgq
        Cgg=self.Cgg
        Cqbq=self.Cqbq
        Cqpq=self.Cqpq
        qf_value=self.qf_value
        gf_value=self.gf_value
        pdf = (self.pdfname,self.setnumber)
        
        match = lambda z,q,f: matching(z,q,f,parton,Cqq,Cqg,Cgq,Cgg,Cqbq,Cqpq,qf_value,gf_value,pdf)
        #match = lambda z,q,f: matchtest(f_value,0,Cgg,z,q,f)
        d = p_map(match,xcopy,qcopy,fcopy)
        self.result = d
        return d
    
    def setxQtable(self,x_value,Q_value):
        self.x_value=x_value
        self.Q_value=Q_value
        
    def setquarkflavourtable(self, qflavour):
        self.qf_value = qflavour
        self.f_value = qflavour
        self.f_value.append(21)
    
    def setpartontype(self,parton):
        # parton should be 'q','g' or 'qg'
        self.parton = parton
    
    def output(self,filename=None):
        if filename is None:
            pass
        else:
            self.filename = filename
            
        self.title1 = " PdfType: replica\n Format: lhagrid1\n ---\n"
        self.title2 = listtotitle(self.x_value)+'\n'
        self.title3 = listtotitle(self.Q_value)+'\n'
        self.title4 = listtotitle(self.f_value)+'\n'
        
        file = open(self.filename, 'w')
        file.write(self.title1)
        file.write(self.title2)
        file.write(self.title3)
        file.write(self.title4)
        for i in range(len(self.x_value)):
            for j in range(len(self.Q_value)):
                for k in range(len(self.f_value)):
                    file.write(str(self.result[k+j*len(self.f_value)+i*len(self.f_value)*len(self.Q_value)]))
                    if k<len(self.f_value):
                        file.write("  ")
                file.write("\n")
            
        file.write(" ---\n")
        file.close()
        
        
        
        
def matching(z,Q,flavour,parton,Cqq,Cqg,Cgq,Cgg,Cqbq,Cqpq,qf_value,gf_value,pdfname):
    pdf = lhapdf.mkPDF(pdfname[0], pdfname[1])
    if parton == 'q' and flavour == 21:
        return 0
    if parton == 'g' and flavour != 21:
        return 0
    if flavour == 21:
        gg = convolute(Cgg,lambda x: PDF(pdf,gf_value,x,Q),z)
        qg = sum([convolute(Cqg,lambda x: PDF(pdf,i,x,Q),z) for i in qf_value])
        return (gg+qg)*z
    else:
        gq = convolute(Cgq, lambda x: PDF(pdf,gf_value,x,Q),z)
        qq = convolute(Cqq, lambda x: PDF(pdf,flavour,x,Q),z)
        qbq = convolute(Cqbq, lambda x: PDF(pdf,-flavour,x,Q),z)
        qpq = 0
        for f in qf_value:
            if f!=flavour and f!=-flavour:
                qpq += convolute(Cqpq, lambda x: PDF(pdf,f,x,Q),z)
        return (gq+qq+qbq+qpq)*z
    
    ### *z is to match the convention of lhapdf
    
def matchtest(t,t2,t3,z,Q,f):
    return z+Q+f


def PDF(pdf,flavour,x,Q):
        return pdf.xfxQ(flavour,x,Q)/x

def convolute(matching_coeffecient, non_perturbative_function,z):
        delta_contribution = matching_coeffecient(1)[0] * non_perturbative_function(1)
        regular_contribution = inte.quad(lambda x:matching_coeffecient(x)[1]*non_perturbative_function(z/x)/x,z,1)[0]
        plus_contribution = 0
        for i in range(2,len(matching_coeffecient(1))):
            n = i-2
            f = lambda x: matching_coeffecient(x)[i] * non_perturbative_function(z/x)/x
            flimit = matching_coeffecient(1)[i] * non_perturbative_function(z)
            integrand = lambda x: (f(x)-flimit)* (np.log(1-x)**(n))/(1-x)
            integral = inte.quad(integrand,z,0.999999999)[0]+f(1)*np.log(1-z)**(n+1)/(n+1) if z<0.999999999 else 0
            plus_contribution += integral
            
        return delta_contribution + regular_contribution + plus_contribution