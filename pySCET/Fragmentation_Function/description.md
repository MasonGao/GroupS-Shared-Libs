# Matching to TMD Fragmentation function from collinear FF

author: Chengtai Tan

## Introduction
The Fragmentation function(or parton distribution function) describes the partonic content of hadrons and is of non-perturbative nature. For example, the quark TMDFF can be defined as

$\mathcal{D}^{\rm bare}_{N_2/q}(z,b_\perp)\equiv\int d^{d-2}k_{2\perp}e^{-ib_\perp\cdot k_{2\perp}\tilde{D}^{\rm bare}_{N_2/q}(z,k_{2\perp})}$

$=\sum_X \frac{1}{z}\int\frac{db_+}{4\pi}e^{ib_+P_{2-}/(2z)}\rm{Tr}<0|\frac{\not{n}}{2}\chi_{\bar{n}}(b_+,0,b_\perp)|N_2(P_2),X><N_2(P_2),X|\bar{\chi}_{\bar{n}(0)}|0>$

For more detailed information, one can refer to arXiv:1908.03831v1(Transverse Parton Distribution and Fragmentation Functions at NNLO: the Quark Case).

For large impact parameter $b_T\sim 1/\Lambda_{\rm QCD}$, FFs are non-perturbative, while in semi-perturbative region, they admits and operator product expansion

$\mathcal{F}_{N/q}^{\rm bare}(z,b_\perp/z,\nu)=z^{2-2\epsilon}\mathcal{D}_{N/q}^{\rm bare}(z,b_\perp,\nu)=\sum_i\int_z^1\frac{d\xi}{\xi}d_{N/i(z/\xi)}\mathcal{C}_{iq}^{\rm bare}(\xi,b_\perp/\xi,\nu)+\mathcal{O}(b_T^2\Lambda_{\rm QCD}^2).$

where $\mathcal{C}_{iq}$ is the perturbative matching coeffecient and $d_{N/i}$ is the collinear FF.

## Numerical Implementation

The scale dependent part of Fragmentation function is encoded in the renormalization group equation. After subtracted by the soft function with $\xi$ prescription. The rapidity depenedence cancels out, the matching coeffecient can be generally written as $\mathcal{C}_{iq}^{(n)}=\sum_{k=0}^nc^{n}_kL_\perp^k$. After convolution, we get, similarly $\mathcal{F}_{iq}^{(n)}=\sum_{k=0}^nf^{n}_kL_\perp^k$. Coeffecients are obtained in programme one by one.


At NLO, the matching coeffecient for TMDFF reads(with $\xi$ prescription implemented)

$\mathcal{C}_{iq}^{(1)}(z,b_\perp/z,L_Q)=\gamma_0^B L_\perp \delta_{iq}\delta(1-z)-P_{iq}^{T(0)}L_\perp+\mathcal{C}_{iq}^{(1)}(z)$

where $\gamma_0^B=3C_F$, $P_{qq}^{(0)}=2C_F[\frac{1+z^2}{(1-z)_+}+\frac{3}{2}\delta(1-z)]$,$P_{gq}^{(0)}=2C_F\frac{1+(1-z)^2}{z}$，$C_{qq}^{(1)}=2C_F(2\frac{1+z^2}{(1-z)_+}\ln(z)+1-z)$,$C_{gq}^{(1)}(z)=2C_F(2\frac{1+(1-z)^2}{z}\ln(z)+z)$


## example code

one can see TMDFF_example.py, TMDFF_minimal for instructions to use the OPEmatching module. 

