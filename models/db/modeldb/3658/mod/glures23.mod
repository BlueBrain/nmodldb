TITLE glutamate summation for estimation of PPR, 2D-3D diffusion
COMMENT
Author: Elena Saftenku, 2003
ENDCOMMENT
NEURON{
POINT_PROCESS GrC_Glures23
RANGE    glu,rPSD,rabs,nu, gluspill,gludir,alpha,Rmf
RANGE taurec,taufacil,tauin, u0, usr,Deff,c0cleft,meandist,rabs,
alpha,nua  
RANGE tm1,td1,ts1
RANGE inclugludir,inclugluspill, Popeak,alpha,Podir,Pospill 
RANGE E1,u1,R1, P1, Nvesicles,PrP, inclugludir1spike,
inclugluspill1spike
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
PARAMETER { Deff=0.052 (um2/ms): effective difusion coefficient 
nu=1.72(1/um2) : density of release sites 
rabs= 3.6 (um) : radius of absorbing boundary
c0cleft = 5.167 (mM): initial glutamate concentration after release
rPSD=0.11 (um): radius of postsynaptic density
taurec=5000 (ms): time constant of recovery, very large in this simulation
tauin=3 (ms): time constant of decay of released substance
taufacil=0(ms): time constant of facilitation
usr=0.3 : vesicle release probability
u0=0 <0,1> : initial value for the "facilitation variable"
meandist=0.24 (um) : minimal limit of spillover [glu] integration
alpha=5 : 1/extracellular volume fraction
h=0.02(um): cleft width
Rmf=2.5(um): radius of mossy fiber
Popeak=0.655 : adjusted open release probability of AMPA receptors
inclugludir=1 : include direct release for the second spike
inclugluspill=1 : include release of spillover glutamate for the second
: spike
inclugludir1spike=1: include direct release for the first spike 
inclugluspill1spike=1: include release of spillover glutamate 
: for the first spike
Nvesicles=2 : number of ready to release vesicles
tm1=0: shift for experimental mEPSC
td1=0: shift for experimental direct EPSC
ts1=0: shift for experimental spillover EPSC
 }
VERBATIM
int i, n,ii, includir[50],incluspill[50];
double l[200000], t0[50],nuac[50],add;
extern float bessj1();
ENDVERBATIM
ASSIGNED{
: Podir
: Pospill 
Nsp : number of spikes
tx1(ms) : time of release
gludir (mM)
gluspill(mM)
vspr
glu (mM)
sum[50] (um)
sum0[50] (um)
sum02[50]
sum2[50]
sum1[50](um2)
sum01[50](um2)
nua (/um2): density of active release sites
R1 : availability of vesicle
u1 : vesicle release probability
E1
P1: release probability
PrP: release probability for the first release
}

INITIAL {
tx1=10000000
glu=0 
gludir=0
gluspill=0
PrP=0
nua=0
R1=1
E1=0
u1=0
Nsp=0
}
BREAKPOINT
{
at_time(tx1)
if (t<=tx1){
glu=0
:Podir=0
:Pospill=0
gludir=0
gluspill=0
}
if (t>tx1){
VERBATIM
l[0]=0;
l[1]=2.4048;l[2]=5.5201;l[3]=8.65;
gluspill=0; gludir=0;
for (ii=1; ii<=n;ii++){
includir[ii]=0;
incluspill[ii]=0;
includir[1]=inclugludir1spike;
includir[2]=inclugludir;
incluspill[1]=inclugluspill1spike;
incluspill[2]=inclugluspill;
sum[ii]=0; i=1; 
 do 
 {if (i>=4) l[i]=PI*(4*i-1)/4;
 sum0[ii]=sum[ii];
 add=(l[i]/rabs)*(l[i]/rabs)*Deff*(t0[ii]-t);
 if (add<-20.0) add=-20.0;
 sum[ii] =sum[ii]+bessj1((l[i]/rabs)*rPSD)/((l[i]/rabs)*
 bessj1(l[i])*
 bessj1(l[i]))* exp(add);
  i++; }
 while (fabs(sum[ii]-sum0[ii])>=0.01);
sum2[ii]=0;i=0;
do
{sum02[ii]=sum2[ii];
sum2[ii]=sum2[ii]+(4/((2*i+1)*PI))*sin((2*i+1)*PI*h/(2*(rabs-Rmf)))*
exp(Deff*(2*i+1)*(2*i+1)*PI*PI*(t0[ii]-t)/(4*(rabs-Rmf)*(rabs-Rmf)));
i++;}
while(fabs(sum2[ii]-sum02[ii])>=0.001);
 gludir =gludir+includir[ii]* 2*c0cleft*rPSD*sum[ii]*sqrt(alpha*sum2[ii])/(rabs*rabs);
if(gludir>c0cleft) gludir=c0cleft;
sum1[ii]=0;i=1;
do
{if (i>=4) l[i]=PI*(4*i-1)/4;
sum01[ii]=sum1[ii];
sum1[ii]=sum1[ii]+(Rmf*bessj1((l[i]/rabs)*Rmf)- meandist 
* bessj1((l[i]/rabs)* meandist))/
((l[i]/rabs)*bessj1(l[i])*bessj1(l[i]))*exp((l[i]/rabs)*(l[i]/rabs)*
Deff*(t0[ii]-t));
i++;}
while (fabs(sum1[ii]-sum01[ii])>=0.0001);
gluspill= gluspill+incluspill[ii]*2*PI*nuac[ii]*c0cleft*rPSD*rPSD*
sum1[ii]*sqrt(sum2[ii]*alpha)/(rabs*rabs);
}
ENDVERBATIM
glu= gludir  + gluspill

: Experimental waveforms
: Podir=(0.94*exp((tx1-t)/0.6(ms))+0.06*exp((tx1-t)/3.57(ms))
: -exp((tx1-t)/0.326(ms)))/0.246*(0.43/0.484)*Popeak
: Pospill=(0.39*exp((tx1-t)/3.25(ms))+0.61*exp((tx1-t)/14.78(ms))-
: exp((tx1-t)/0.721(ms)))/0.682*(0.125/0.484)*Popeak
}
}
NET_RECEIVE (weight,Eav,R, u,P,nspike, tsyn (ms))
{
INITIAL
{
R=1
Eav=0
u=u0
tsyn=t
nspike=0
}
nspike=nspike+1
vspr=((1-R-Eav)/taurec+(R-1)/tauin)/(1/tauin-1/taurec)
   R=1+exp((tsyn-t)/taurec)*vspr+exp((tsyn-t)/tauin)*(R-1-vspr)
   Eav=Eav*exp((tsyn-t)/tauin)
  if (taufacil>0){
   u=u*exp((tsyn-t)/taufacil)
   }else {
   u=usr
    }
if (taufacil>0) {
state_discontinuity (u, u+usr*(1-u))  
}
if (Nvesicles==1){P=u}
if(Nvesicles==2){P=(1-(1-u)^2)}
nua=0
: printf("Eav= %f\n",Eav)
if (nspike==1 || Nvesicles==1){
state_discontinuity (Eav, Eav+P*R)
state_discontinuity (nua, nua+nu*P*R)
}
if(nspike==2&& Nvesicles==2){
state_discontinuity (Eav, Eav+u*PrP+P*R)
state_discontinuity(nua, nua+nu*(u*PrP+P*R))
}
:printf("R= %f\n",R)
:printf("u= %f\n",u)
{state_discontinuity(R, R-P*R)}
PrP=P
tsyn=t
tx1=t 
E1=Eav
R1=R
u1=u
P1=P
Nsp=nspike
VERBATIM
n=Nsp;
t0[n]=tx1;
nuac[n]=nua;
ENDVERBATIM
}

