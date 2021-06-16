
COMMENT

This file is based on na16.mod of Hu et al. (nn, 2009) 

******************************************************
Sodium channel, Hodgkin-Huxley style kinetics.  

Kinetics were fit to data from Huguenard et al. (1988) and Hamill et al. (1991)

qi is not well constrained by the data, since there are no points between -80 and -55.  So this was fixed at 5 
while the thi1, thi2, Rg, Rd were optimized using a simplex least square proc

voltage dependencies are shifted approximately from the best fit to give higher threshold

Author: Zach Mainen, Salk Institute, 1994, zach@salk.edu
******************************************************

SHL modified as follows:
 1-1) qa: 6 -> 8
 1-2) qinf: 6.2 -> 8 mV, since qa and qinf = 8 mV  (Kim Jonas nn 12).
 cf. qa and qinf = 5.8 and 7.1 in Hu Shu nn 09
 2) tha, -43 -> -38 (Kim Jonas nn 12)
 3) use cnexp as an integration mtd
SHL End

ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX Na
	USEION na READ ena WRITE ina
	RANGE m, h, gna, gbar
	GLOBAL tha, thi1, thi2, qa, qi, qinf, thinf
	RANGE minf, hinf, mtau, htau
	GLOBAL Ra, Rb, Rd, Rg
	GLOBAL q10, temp, tadj, vmin, vmax, vshift
}

PARAMETER {
	gbar = 0.01   	(mS/cm2)
	vshift = 0	(mV)		: voltage shift (affects all)

	tha  = -38	(mV)		: v 1/2 for act curve
	qa   = 8	(mV)		: act slope		
	Ra   = 0.182	(/ms)	: prop const of alpha_m (/ms)	
	Rb   = 0.124	(/ms)	: prop const of beta_m (/ms)	

	thi1  = -50	(mV)		: v shift of alpha_h	
	thi2  = -75	(mV)		: v shift of beta_h	
	qi   = 5	(mV)	    : inact tau slope
	thinf  = -73	(mV)	: v1/2 of inact curve	
	qinf  = 8	(mV)		: Slope of inact curve; 6.2 -> 8
	Rg   = 0.0091	(/ms)	: prop const of beta_h(v)	
	Rd   = 0.024	(/ms)	: prop const of alpha_h(v)

	temp = 30	(degC)		: original temp 
	q10  = 2.3			: temperature sensitivity

	v 		(mV)
	dt		(ms)
	celsius		(degC)
	vmin = -120	(mV)
	vmax = 100	(mV)
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

ASSIGNED {
	ina 	(mA/cm2)
	gna		(pS/um2)
	ena		(mV)
	minf 		hinf
	mtau (ms)	htau (ms)
	tadj
}
 

STATE { m h }

INITIAL { 
	rates(v+vshift)
	m = minf
	h = hinf
	tadj = q10^((celsius - temp)/10)
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    gna = tadj*gbar*m*m*m*h
	ina = gna * (v - ena)
} 

DERIVATIVE states {
    rates(v+vshift)      :   at the current v and dt.   
    m' = (minf-m)*tadj/mtau
    h' = (hinf-h)*tadj/htau
}



PROCEDURE rates(vm) {  
    LOCAL  a, b

	a = trap0(vm,tha,Ra,qa)
	b = trap0(-vm,-tha,Rb,qa)
	mtau = 1/(a+b)
	minf = a*mtau

	:"h" inactivation 

	a = trap0(vm,thi1,Rd,qi)
	b = trap0(-vm,-thi2,Rg,qi)
	htau = 1/(a+b)
	hinf = 1/(1+exp((vm-thinf)/qinf))
}


FUNCTION trap0(v,th,a,q) {
	if (fabs(v-th) > 1e-6) {
	    trap0 = a * (v - th) / (1 - exp(-(v - th)/q))
	} else {
	    trap0 = a * q
 	}
}	
