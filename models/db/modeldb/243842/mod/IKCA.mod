
TITLE Calcium activated Potassium channel (IK)
  
: Author: Chitaranjan Mahapatra (chitaranjan@iitb.ac.in)
: Computational Neurophysiology Lab
: Indian Institute of Technology Bombay, India 

: For details refer: 
: Mahapatra C, Brain KL, Manchanda R, A biophysically constrained computational model of the action potential 
: of mouse urinary bladder smooth muscle. PLOS One (2018) 


NEURON {
    SUFFIX IKCA
    USEION k READ ek WRITE ik
    USEION ca READ cai
    RANGE gkbar, gk,ik,ph,pm,vah,cvh,ckh,mt,ht
	RANGE vhalf ,cal,tau,minf,taum,linf,taul,slp,min, q10
	 
	       }

UNITS {
	(molar) = (1/liter)
	(mM)	= (millimolar)
	(S)  	= (siemens)
	(mA) 	= (milliamp)
	(mV) 	= (millivolt)
}

PARAMETER {
	
	gkbar = 0.007 (S/cm2) 
	ca0 =  0.00015 	(mM)   
	cal = 0 (mM)
	tau = 2		(ms)	   
	mt = 1
	ht = 1
	vha = -38
	slp = 15
   	
    pm = 2
	ph = 1
	:ap  = 1.5
	cvh = -66 (mV)
	ckh = -8 (mV)
	min = 0.27
	
	celsius = 37(degC)
	q10 = 2.3
		
	}


ASSIGNED {
	ik (mA/cm2)
	gk (S/cm2)   
	minf
	taum (ms)
	linf (mv)
	taul (ms)
	vhalf (mV)
	v       (mV)
	cai  (mM)
	ek		(mV)
}


STATE {
	m   
	l 
	ca_i (mM)
}

BREAKPOINT {
	SOLVE states METHOD cnexp
     gk = gkbar * l^pm * m^ph	 
	 ik = gk * (v - ek)
	 cal = ca_i
	 }
	 
DERIVATIVE states {
	ca_i' = -(ca_i-ca0)/tau
	rates(v,ca_i)
	m' = (minf-m)/taum	
	l' = (linf-l)/taul
	
}

INITIAL {
    ca_i = ca0
	rates(v, ca_i)
	m = minf
	l = linf
}

PROCEDURE rates(v (mV), c (mM)) {
	LOCAL a,qt
        qt=q10^((celsius-37)/10)
	linf =  vha-152 + (59.2*exp(-.09*c*1e3)) + (96.7*exp(-.47*c*1e3))
    linf = min/(1+exp((linf-v)/(slp)))
	
	
	taul = mt*5*(17/(1+(((v+20.52)/35)^2)))
	
	minf = 0.32/(1 + exp((v-cvh)/ckh))
	taum = ht * qt* ( 0.5*(90.9699*(1.0-(1.0/(1.0+exp((v+13.9629)/45.3782)))*(1.0/(1.0+exp(-(v+9.49866)/3.3945)))))	)
	
}

FUNCTION trap0(v,th,a,q) {
	if (fabs(v-th) > 1e-6) {
	       trap0 = a * (v - th) / (1 - exp(-(v - th)/q))
	} else {
	        trap0 = a * q
 	}
}	