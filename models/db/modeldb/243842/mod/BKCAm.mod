
TITLE Calcium activated Potassium channel (BK)  

: Author: Chitaranjan Mahapatra (chitaranjan@iitb.ac.in)
: Computational Neurophysiology Lab
: Indian Institute of Technology Bombay, India 

: For details refer: 
: Mahapatra C, Brain KL, Manchanda R, A biophysically constrained computational model of the action potential 
: of mouse urinary bladder smooth muscle. PLOS One (2018) 


NEURON{
SUFFIX BKCAm
USEION k READ ek WRITE ik
USEION ca READ cai 
 RANGE  ik, gkbar, kon,kcoff, kooff, c0o0c, c0c1c, o0c0c,o0o1c, o4c4c,c4o4c, cai, a, b, hva , hvi,sla,sli,o
}

UNITS {
	(molar) = (1/liter)
	(mM)	= (millimolar)
	(S)  	= (siemens)
	(mA) 	= (milliamp)
	(mV) 	= (millivolt)
}

PARAMETER{
gkbar = 0.024 (S/cm2)
kon = 40633
kcoff = 11
kooff = 1.1

c0o0c = 0.02162
c0c1c = 4

o0c0c = 318.1084
o0o1c = 4

o4c4c = 0.35
c4o4c = 0.002

hva = 35
hvi = 15.776
sla = 380
sli = 330


cai	= 0.0001(mM)
}

ASSIGNED{

v (mV)
ek 	(mV)
ik	(mA/cm2)
a
b
o
	

c0c1 (/ms)
c1c0 (/ms)
c1c2 (/ms)
c2c1 (/ms)
c2c3 (/ms)
c3c2 (/ms)
c3c4 (/ms)
c4c3 (/ms)

o0o1 (/ms)
o1o0 (/ms)
o1o2 (/ms)
o2o1 (/ms)
o2o3 (/ms)
o3o2 (/ms)
o3o4 (/ms)
o4o3 (/ms)

c0o0 (/ms)
o0c0 (/ms)

c1o1 (/ms)
o1c1 (/ms)

c2o2 (/ms)
o2c2 (/ms)

c3o3 (/ms)
o3c3 (/ms)

c4o4 (/ms)
o4c4 (/ms)

}


STATE {c0 c1 c2 c3 c4 o0 o1 o2 o3 o4 }

BREAKPOINT {SOLVE kin METHOD sparse
o = o0 + o1 + o2 + o3 + o4



ik = gkbar * o * ( v - ek )
}

 INITIAL { 
 
SOLVE kin STEADYSTATE sparse

}

KINETIC kin {

       rates(v)

~ c0<->c1 (c0c1, c1c0)
~ c1<->c2 (c1c2, c2c1)
~ c2<->c3 (c2c3, c3c2)
~ c3<->c4 (c3c4, c4c3)

~ o0<->o1 (o0o1, o1o0)
~ o1<->o2 (o1o2, o2o1)
~ o2<->o3 (o2o3, o3o2)
~ o3<->o4 (o3o4, o4o3)

~ c0<->o0 (c0o0, o0c0)
~ c1<->o1 (c1o1, o1c1)
~ c2<->o2 (c2o2, o2c2)
~ c3<->o3 (c3o3, o3c3)
~ c4<->o4 (c4o4, o4c4)

CONSERVE c0 + c1 + c2 + c3 + c4 + o0 + o1 + o2 + o3 + o4 = 1
}

PROCEDURE rates(v(mV)){
UNITSOFF
a = exp (hva*v/sla)
b = exp (hvi*v/sli)

c0o0 = c0o0c   * a
c1o1 = 0.000869  * a
c2o2 = 0.0000281 * a
c3o3 = 0.000781  * a
c4o4 = c4o4c * a

o0c0 = o0c0c *  b
o1c1 = 144.1736 * b
o2c2 = 32.6594 * b
o3c3 = 0.095312 * b
o4c4 = o4c4c * b
UNITSON
}

PROCEDURE prates(cai(mM)){

UNITSOFF

c0c1 = c0c1c * kon * cai
c1c2 = 3 * kon * cai
c2c3 = 2 * kon * cai
c3c4 = kon * cai

c4c3 = 4 * kcoff * cai
c3c2 = 3 * kcoff * cai
c2c1 = 2 * kcoff * cai
c1c0 =  kcoff * cai

o0o1 = 4 * kon * cai
o1o2 = 3 * kon * cai
o2o3 = 2 * kon * cai
o3o4 =  kon * cai

o4o3 = 4 * kooff * cai
o3o2 = 3 * kooff * cai
o2o1 = 2 * kooff * cai
o1o0 =  kooff * cai

UNITSON
}