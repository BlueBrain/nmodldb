TITLE Long-term potentiation and depression
COMMENT
Realization of Migliore and Lansky model of long-term potentiation and depression (1999). In contrast to the model of Migliore and Lansky, an activation of autocatalytic processes is controlled not by
postsynaptic depolarization, but by [Ca] influx through the NMDA 
receptor channels. 
Author: Elena Saftenku, 2001
ENDCOMMENT
NEURON{
POINT_PROCESS  Gr_LTP1
POINTER picanmda
RANGE gamma,eta,nu1,nu2,pp,pd,gdel1,gdel2,Mp,Md,Ap,Ad
}
UNITS{
(nA)=(nanoamp)
(mA)=(milliamp)
(mV)=(millivolt)
}
PARAMETER {
 gamma=0.34 (/ms) : rate constant of the production of 
:[Ca]-dependent protein
 eta=0.003 (/ms) : rate constant of protein degradation 
 nu1=0.065 (/ms) : rate constant of activation of autocatalytic 
: process Np
 nu2=0.065 (/ms): the same for Nd
:Parameters of autocatalytic processes
 pp=0.0000095 (/ms)
 pd=0.00019(/ms)
 gdel1= 2.4(/nA-ms)
 gdel2= 2.4(/nA-ms)
 Mp= 3e-5(nA/ms)
 Md=3e-4(nA/ms)
 Ap=1.625(nA2)
 Ad=0.55(nA2)
}
ASSIGNED{
picanmda (nA): reference to Ca2+ current through NMDA receptors
 }
STATE{
Np (nA) 
Nd (nA)  
messenger (nA) <1e-04>
}
INITIAL {
messenger=0
  Nd=0
  Np=0
}
BREAKPOINT
{
SOLVE states METHOD derivimplicit
}

DERIVATIVE states {
messenger' = -gamma*picanmda - eta*messenger
  Np' = nu1*messenger  - (pp - picanmda*gdel1)*Np +(Mp*Np*Np)/(Ap+Np*Np)
  Nd' = nu2*messenger  -(pd - picanmda*gdel2)*Nd  +(Md*Nd*Nd)/(Ad+Nd*Nd)
}

