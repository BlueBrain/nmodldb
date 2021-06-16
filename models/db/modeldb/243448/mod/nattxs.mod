TITLE TTX Sensitive Current for bladder small DRG neuron soma model
: Author: Darshan Mandge (darshanmandge@iitb.ac.in)
: Computational Neurophysiology Lab
: Indian Institute of Technology Bombay, India 

: For details refer: 
: A biophysically detailed computational model of bladder small DRG neuron soma 
: Darshan Mandge and Rohit Manchanda, PLOS Computational Biology (2018)

NEURON {
	SUFFIX nattxs
	USEION na READ ena WRITE ina
	RANGE gbar, ena, ina
	RANGE minf, hinf, mtau, htau
	THREADSAFE
}

UNITS {
	(S)  = (siemens)
	(mV) = (millivolts)
	(mA) = (milliamp)
}

PARAMETER {
	gbar = 0.0001 (S/cm2)
}

ASSIGNED {
	v		(mV) 
	ena 	(mV)
	ina		(mA/cm2)
	g		(S/cm2)
	
	htau	(ms)
	mtau	(ms)
	minf
	hinf
}

STATE { m h }

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gbar*m*m*m*h
	ina = g * (v-ena)
}

INITIAL {
	rates(v)
	m = minf
	h = hinf
}

DERIVATIVE states {
	rates(v)
	m' = (minf - m)/mtau
	h' = (hinf - h)/htau
}

PROCEDURE rates(v (mV)) {
	LOCAL alpha_m, beta_m, alpha_h, beta_h
	
	UNITSOFF
	minf 	= 1/(1+exp(((-25.8)-v)/7.8)) 		: Data Fit: Yoshimura et al., 1996 
	hinf 	= 1/(1+exp((v+55.8)/8.9))			: Data Fit: Yoshimura et al., 1996
	
	
	alpha_m = 15.5/(1+exp((v-5)/(-12.08)))
	beta_m  = 35.2/(1+exp((v+72.7)/16.7)) 
	mtau 	= 1 / (alpha_m+beta_m)				: Adapted from Sheets et al.,2007 Nav1p7 Channel
	
	
	alpha_h = (0.23688)*exp(-(v+115)/46.33)
	beta_h 	= (10.8/2.5)/(1+exp((v-11.8)/-11.998))
	htau 	= 1 / (alpha_h + beta_h) 			: Adapted From Baker, 2005 nattxs.mod Channel

	UNITSON
}