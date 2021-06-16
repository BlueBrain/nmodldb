TITLE 3D diffusion, closed boundary
COMMENT
Author: Elena Saftenku, 2003
ENDCOMMENT
NEURON{
POINT_PROCESS GrC_Glubes4
RANGE    glu,rPSD,rabs,nu, gluspill,gludir,alpha,Rmf
RANGE Deff,c0cleft,meandist,rabs,alpha   
RANGE inclugludir,inclugluspill, Popeak,alpha,Podir,Pospill
RANGE td1,tm1,ts1
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
PARAMETER { Deff=0.043 (um2/ms) :effective diffusion coefficient
nu=2 (/um2):density of release sites  
rabs= 4.4 (um) :radius of absorbing boundary
c0cleft = 8.769 (mM):initial [glu]in the cleft after vesicle release
rPSD=0.11 (um): radius of postsynaptic density
meandist=0.29 (um): the lower limit of spillover [glu] integration
Rmf=2.9(um): radius of mossy fiber terminal
alpha=5: 1/extracellular volume fraction
h=0.02 (um):synaptic cleft width
Popeak=0.678: adjusted peak open probability of AMPA receptors
inclugludir=1 : inclusion of direct component
inclugluspill=1: inclusion of spillover component
td1=0 (ms) : 0.09 (ms) , shift of experimental mEPSC
ts1=0 (ms) : 0.16 (ms) , shift of experimental direct EPSC
tm1=0 (ms) : 0.16 (ms) , shift of experimental spillover EPSC
}
VERBATIM
int i;
double l[100000];
float bessj1();
float bessj0();
ENDVERBATIM
ASSIGNED{
Podir
Pospill 
tx1(ms)
gludir (mM)
gluspill(mM)
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
}
if(t>tx1) {
VERBATIM
l[0]=0;
l[1]=3.8317;l[2]=7.01558667;l[3]=10.1737;
sum=0; i=1; 
do 
{if (i>=4) {l[i]=0;
l[i]=PI*(4*i+1)/4;
}
sum0=sum;
sum =sum+bessj1((l[i]/rabs)*rPSD)/((l[i]/rabs)*bessj0(l[i])*bessj0(l[i]))
* exp((l[i]/rabs)*(l[i]/rabs)*Deff*(tx1-t));
 i++; }
while (fabs(sum-sum0)>=0.01);
sum2=0;i=1;
do
{sum02=sum2;
sum2=sum2+(2/(i*PI))*sin(i*PI*h/(rabs-Rmf))*
exp(Deff*i*i*PI*PI*(tx1-t)/((rabs-Rmf)*(rabs-Rmf)));
i++;}
while(fabs(sum2-sum02)>=0.00001);
ENDVERBATIM
gludir= 2*c0cleft*rPSD*((rPSD/2)+sum)*alpha*
((h/(rabs-Rmf))+sum2)/(rabs*rabs)
if(gludir>c0cleft){gludir=c0cleft}
VERBATIM
sum1=0;i=1;
do
{if (i>=4) { 
l[i]=0;
l[i]=PI*(4*i+1)/4;
}
sum01=sum1;
sum1= sum1+(Rmf*bessj1((l[i]/rabs)*Rmf)- meandist 
*bessj1((l[i]/rabs)*meandist))/
((l[i]/rabs)*bessj0(l[i])*bessj0(l[i]))*exp((l[i]/rabs)*(l[i]/rabs)
*Deff*(tx1-t));
i++;}
while(fabs(sum1-sum01)>=0.0001);
ENDVERBATIM
gluspill=2*PI*nu*c0cleft*rPSD*rPSD*((Rmf*Rmf-meandist*meandist)/2+sum1)*((h/(rabs-Rmf))+sum2)*alpha/(rabs*rabs)
glu= inclugludir*gludir  +inclugluspill*gluspill

:Experimental waveforms
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

VERBATIM
float bessj1(x)
float x;
{
float ax,z;
double xx,y,ans,ans1,ans2;
if ((ax=fabs(x))<8.0) {
y=x*x; 
ans1=x*(72362614232.0+y*(-7895059235.0+y*(242396853.1+
  y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
  ans2=144725228442.0+y*(2300535178.0+y*(18583304.74+
  y*(99447.43394+y*(376.9991397+y*1.0))));
  ans=ans1/ans2; 
  } else {
  z=8.0/ax; y=z*z; xx=ax-2.356194491;
  ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4+y*(0.2457520174e-5+
  y*(-0.240337019e-6))));
  ans2=0.04687499995+y*(-0.2002690873e-3+y*(0.8449199096e-5
  +y*(-0.88228987e-6+y*0.105787412e-6)));
  ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
if (x<0.0) ans=-ans;
}
return ans;
}
float bessj0(x)
float x;
{
float ax,z;
double xx,y,ans,ans1,ans2;
if ((ax=fabs(x))<8.0) {
y=x*x; 
ans1=57568490574.0+y*(-13362590354.0+y*(651619640.7+
  y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
  ans2=57568490411.0+y*(1029532985.0+y*(9494680.718+
  y*(59272.64853+y*(267.8532712+y*1.0)))); ans=ans1/ans2;
  } else {
  z=8.0/ax;y=z*z;  xx=ax-0.785398164;
  ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
  +y*(-0.2073370639e-5+y*0.2093887211e-6)));
  ans2=-0.1562499995e-1+y*(0.1430488765e-3+
  y*(-0.6911147651e-5+y*(0.7621095161e-6-y*0.934945152e-7)));
  ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
  }
return ans;
}
ENDVERBATIM

