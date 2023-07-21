from numpy import pi
from scipy.special import zeta
import numpy as np

CA=3.
CF=4./3.
TF=0.5
nf=5.
Nc=3


# cusp anomalous dimension
# convention: gamma_cusp^0=4, Gamma_cusp^0=4CF
gamma_cusp=np.array([4., \
4*((67/9-pi**2/3)*CA-20/9*TF*nf), \
4*(CA**2 *(245/6 - 134/27*(pi**2) +11/45*(pi**4) + 22/3*zeta(3))+CA*TF*nf*(-418/27 + 40/27*(pi**2) - 56/3*zeta(3)) + CF*TF*nf*(-55/3 + 16*zeta(3)) - 16/27*(TF*nf)**2), \
1553/CF])

Gamma_cusp = gamma_cusp * CF

# non-cusp anomalous dimension for Hard function
gamma_mu=np.array([-6*CF, \
CF**2*(-3+4*(pi**2)-48*zeta(3))  +  CF*CA*(-961/27-11/3*(pi**2)+52*zeta(3))  +  CF*TF*nf*(260/27+4/3*(pi**2)),   \
CF**3 * (-29 - 6*pi**2 - 16 *pi**4/5 - 136 * zeta(3) + 32*pi**2/3*zeta(3) + 480 * zeta(5)) + CF**2 * CA * ( - 151/2 + 410*pi**2/9 + 494*pi**4/135 - 1688*zeta(3)/3 - 16*pi**2*zeta(3)/3 - 240 *zeta(5)) + CF * CA**2 * (-139345/1458 - 7163*pi**2/243 - 83*pi**4/45 + 7052*zeta(3)/9 - 88*pi**2*zeta(3)/9 - 272*zeta(5)) + CF**2 * TF * nf * (5906/27 - 52*pi**2/9 - 56*pi**4/27 + 1024*zeta(3)/9) + CF * CA * TF * nf * (- 34636/729 + 5188*pi**2/243 + 44*pi**4/45 - 3856*zeta(3)/27) + CF * TF**2 * nf**2* (19336/729 - 80*pi**2/27 - 64/27*zeta(3)) ] )


# QCD beta function
beta=np.array([11/3*CA-4/3*TF*nf, \
    34/3*CA**2-20/3*CA*TF*nf-4*CF*TF*nf, \
    2857/54*(CA**3) + (2*(CF**2) - 205/9*CF*CA - 1415/27*(CA**2))*TF*nf + (44/9*CF + 158/27*CA)*(TF*nf)**2, \
    149753/6 + 3564*zeta(3) - nf*(1078361/162 + 6508/27 *zeta(3)) + nf**2 *(50065/162 + 6472*zeta(3)/81) + nf**3 *1093/729])
