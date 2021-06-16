TITLE 3D diffusion with absorbing boundary
COMMENT
Author: Elena Saftenku, 2003
ENDCOMMENT
NEURON{
POINT_PROCESS GrC_Glubes3
RANGE glu,rPSD,rabs,nu, gluspill,gludir,alpha,Rmf
RANGE Deff,meandist,rabs,alpha  
RANGE inclugludir,inclugluspill, Popeak,alpha,Podir,Pospill, Podir1, Posum, Posum1,td1,tm1,ts1
}
UNITS{
(molar)=(1/liter)
(mM)=(millimolar)
(um)=(micron)
(nA)=(nanoamp)
}
CONSTANT {
PI=3.1415927
}
PARAMETER { Deff=0.043 (um2/ms) : effective diffusion coefficient
nu=2 (/um2)  :   density of release sites
rabs= 4.4 (um) : radius of absorbing boundary
c0cleft = 8.769 (mM): initial[glu] in the cleft after release
rPSD=0.11 (um): radius of postsynaptic density
meandist=0.29 (um) : minimal limit of spillover glutamate integration
Rmf=2.9(um): radius of mossy fiber
alpha=5: 1/extracellular volume fraction
h=0.02 (um): synaptic cleft width
Popeak=0.678 : adjusted peak open probability of AMPA receptors
inclugludir=1 : inclusion of direct component
inclugluspill=1 : inclusion of spillover component
tm1=0 (ms) : 0.09 (ms) , shift of experimental mEPSC
td1=0 (ms) : 0.16(ms), shift of experimental direct EPSC 
ts1=0 (ms): 0.16 (ms), shift of experimental spillover EPSC 
}
VERBATIM
int i;
double l[100000];
extern float bessj1();
ENDVERBATIM
ASSIGNED{
Podir
Pospill 
Podir1
Posum
Posum1
tx1(ms)
gludir (mM)
gluspill(mM)
vspr
glu (mM)
sum (um)
sum0 (um)
sum02
sum2
sum1(um2)
sum01(um2)
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
Posum=0
Posum1=0
Podir1=0
}
if(t>tx1) {
VERBATIM
l[0]=0;
l[1]=2.4048;l[2]=5.5201;l[3]=8.6535;
sum=0; i=1; 
do 
{if (i>=4) l[i]=PI*(4*i-1)/4;
sum0=sum;
sum =sum+bessj1((l[i]/rabs)*rPSD)/((l[i]/rabs)*bessj1(l[i])*bessj1(l[i]))
* exp((l[i]/rabs)*(l[i]/rabs)*Deff*(tx1-t));
 i++; }
while (fabs(sum-sum0)>=0.01);
sum2=0;i=0;
do
{sum02=sum2;
sum2=sum2+(4/((2*i+1)*PI))*sin((2*i+1)*PI*h/(2*(rabs-Rmf)))*
exp(Deff*(2*i+1)*(2*i+1)*PI*PI*(tx1-t)/(4*(rabs-Rmf)*(rabs-Rmf)));
i++;}
while(fabs(sum2-sum02)>=0.001);
ENDVERBATIM
gludir= 2*c0cleft*rPSD*sum*alpha*sum2/(rabs*rabs)
if (gludir>c0cleft){gludir=c0cleft}
VERBATIM
sum1=0;i=1;
do
{if (i>=4) l[i]=PI*(4*i-1)/4;
sum01=sum1;
sum1= sum1+(Rmf*bessj1((l[i]/rabs)*Rmf)- meandist 
*bessj1((l[i]/rabs)*meandist))/
((l[i]/rabs)*bessj1(l[i])*bessj1(l[i]))*exp((l[i]/rabs)*(l[i]/rabs)
*Deff*(tx1-t));
i++;}
while(fabs(sum1-sum01)>=0.001);
ENDVERBATIM
gluspill=2*PI*nu*c0cleft*rPSD*rPSD*sum1*sum2*alpha/(rabs*rabs)
: printf("gluspill= %f\n", gluspill)
glu= inclugludir*gludir  +inclugluspill*gluspill

: Experimental waveforms
Podir=(0.94*exp((tx1-td1-t)/0.37(ms))+0.06*exp((tx1-td1-t)/2.2(ms))
  -exp((tx1-td1-t)/0.199(ms)))/0.249*(0.43/0.484)*Popeak :direct EPSC
Pospill=(0.39*exp((tx1-ts1-t)/2.0(ms))+0.61*exp((tx1-ts1-t)/9.1(ms))-
 exp((tx1-ts1-t)/0.44(ms)))/0.682*(0.125/0.484)*Popeak: spillover EPSC
Podir1=(0.94*exp((tx1-tm1-t)/0.3(ms))+0.06*exp((tx1-tm1-t)/3.1(ms))
  -exp((tx1-tm1-t)/0.12(ms)))/0.35*Popeak : mEPSC
Posum= Podir + Pospill
Posum1=Podir1 + Pospill
}
}
NET_RECEIVE (weight)
{
tx1=t 
}

