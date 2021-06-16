TITLE 2D-3D diffusion without boundary
COMMENT
Author: Elena Saftenku, 2003
ENDCOMMENT
NEURON{
POINT_PROCESS GrC_Gludif23
RANGE    PRE, glu,rPSD,h,nu,gludir,gluspill
RANGE Deff,meandist,rabs,h,Rmf  
RANGE inclugludir,inclugluspill, Popeak,alpha,Podir,Pospill
RANGE td1,ts1,tm1
}
UNITS{
(molar)=(1/liter)
(mM)=(millimolar)
(um)=(micron)
(nA)=(nanoamp)
PI=(pi) (1)
}
PARAMETER {
nu=0.94(/um2) :density of release sites
rabs=0(um) :radius of absorbing boundary
h=0.02 (um): cleft width
alpha=5 : 1/exstracellular volume fraction
Deff=0.08 (um2/ms): effective diffusion coefficient of glutamate
c0cleft = 8.769 (mM): initial glutamate concentration in the cleft
rPSD=0.11 (um): radius of postsynaptic density
meandist=0.29 (um)
Rmf=2.9(um): radius of mossy fiber terminal
Popeak=0.662: adjusted peak open probability of AMPA receptors
inclugludir=1: inclusion of direct component
inclugluspill=1 : inclusion of spillover component 
td1=0 (ms) : 0.09 (ms) , shift of experimental mEPSC
ts1=0 (ms) : 0.16(ms) , shift of direct EPSC
tm1= 0 (ms) : 0.16 (ms) , shift of spillover EPSC
 }

ASSIGNED{
Podir
Pospill
tx1(ms)
gludir (mM)
gluspill(mM)
glu (mM)
}
INITIAL {
tx1=10000000
glu=0 
gludir=0
gluspill=0
}

BREAKPOINT
{
at_time(tx1)
if (t<=tx1){
glu=0
gludir=0
gluspill=0
Podir=0
Pospill=0
}
if(t>tx1) {
UNITSOFF
gludir= sqrt(2)*c0cleft*sqrt(h)*sqrt(alpha)/sqrt(sqrt(4*PI*Deff*
(t-tx1)))*(1-exp(rPSD*rPSD/(4*Deff*(tx1-t))))
if (gludir>c0cleft){gludir=c0cleft}
gluspill = sqrt(2)*nu*c0cleft*sqrt(h)*rPSD*rPSD*PI*sqrt(alpha)*
(1/sqrt(sqrt(4*PI*Deff*(t-tx1))))*(exp(meandist*meandist/
(4*Deff*(tx1-t)))-exp(Rmf*Rmf/(4*Deff*(tx1-t))))
UNITSON
glu= inclugludir*gludir  +inclugluspill*gluspill

: Experimental waveforms
Podir=(0.94*exp((tx1-t)/0.37(ms))+0.06*exp((tx1-t)/2.2(ms))
-exp((tx1-t)/0.199(ms)))/0.249*(0.43/0.484)*Popeak
Pospill=(0.39*exp((tx1-t)/2.0(ms))+0.61*exp((tx1-t)/9.1(ms))-
exp((tx1-t)/0.44(ms)))/0.682*(0.125/0.484)*Popeak
}
}
NET_RECEIVE (weight)
{
tx1=t 
}

