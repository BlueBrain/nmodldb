TITLE 2D diffusion without boundary
COMMENT
Author: Elena SAFTENKU, 2003
ENDCOMMENT
NEURON{
POINT_PROCESS GrC_Gludif2
RANGE    PRE, glu,rPSD,nu,gludir,gluspill
RANGE Deff,meandist,rabs,Rmf 
RANGE inclugludir,inclugluspill, Popeak,alpha,Podir,Pospill, Podir1
RANGE td1,tm1,ts1
}

UNITS{
(molar)=(1/liter)
(mM)=(millimolar)
(um)=(micron)
(nA)=(nanoamp)
PI=(pi)  (1)
}

PARAMETER {
alpha=1   :extracellular volume fraction
nu=0.33(/um2)  :density of release sites
rabs=0 (um)     :absorbing bpundary
Deff=0.2 (um2/ms):efffective glutamate concentration
c0cleft = 8.769 (mM) :initial glutamate concentration
rPSD=0.11 (um)   :radius of PSD
meandist=0.23 (um): lowest limit of integration
Rmf=2.8 (um) :radius of mossy fiber terminal
Popeak=0.634 : adjusted peak open probability of AMPA receptors
inclugludir=1: inclusion of direct component
inclugluspill=1: inclusion of spillover component
td1=0 (ms) : 0.16(ms) , shift of direct EPSC
tm1=0 (ms) : 0.09 (ms) , shift of experimental mEPSC
ts1=0 (ms) : 0.16 (ms) , shift of spillover EPSC
}
ASSIGNED{
Podir
Podir1
Pospill 
tx1(ms)
gludir (mM)
gluspill(mM)
vspr
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
Podir1=0
Pospill=0
}
if(t>tx1) {
gludir= c0cleft*(1-exp(rPSD*rPSD/(4*Deff*(tx1-t))))
if(gludir>c0cleft){gludir=c0cleft}
gluspill=PI*nu*c0cleft*rPSD*rPSD*
(exp(meandist*meandist/(4*Deff*(tx1-t)))-exp(Rmf*Rmf/(4*Deff*(tx1-t))))
 glu= inclugludir*gludir  +inclugluspill*gluspill
: glu = (c0cleft*rPSD*rPSD/(4*Deff*(t-tx1)))*
: exp((0.46(um)*0.46(um))/(4*Deff*(tx1-t)))

: Experimental waveforms
Podir=(0.94*exp((tx1-t)/0.37(ms))+0.06*exp((tx1-t)/2.2(ms))
  -exp((tx1-t)/0.199(ms)))/0.249*(0.43/0.484)*Popeak
Podir1=(0.94*exp((tx1-t)/0.3(ms))+0.06*exp((tx1-t)/3.1(ms))
  -exp((tx1-t)/0.12(ms)))/0.35*Popeak :           mEPSC
Pospill=(0.39*exp((tx1-t)/2.0(ms))+0.61*exp((tx1-t)/9.1(ms))-
 exp((tx1-t)/0.44(ms)))/0.682*(0.125/0.484)*Popeak
}
}
NET_RECEIVE (weight)
{
tx1=t 
}

