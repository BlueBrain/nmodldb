TITLE BKCa current for bladder small DRG neuron soma model
: Author: Darshan Mandge (darshanmandge@iitb.ac.in)
: Computational Neurophysiology Lab
: Indian Institute of Technology Bombay, India 

: For details refer: 
: A biophysically detailed computational model of bladder small DRG neuron soma 
: Darshan Mandge and Rohit Manchanda, PLOS Computational Biology (2018)

NEURON {
	SUFFIX bkca
	USEION k READ ek WRITE ik
	USEION ca READ cai
	RANGE g, ik
	RANGE gbar, ninf, ntau, vhalf, sf1, pCa
}

UNITS {
	(molar) = (1/liter)
	(mM)	= (millimolar)
	(S)  	= (siemens)
	(mA) 	= (milliamp)
	(mV) 	= (millivolt)
}

PARAMETER {
        gbar	= 0.0009 (S/cm2)
		ek 				 (mV)
}

ASSIGNED {
        v       (mV)
        ik		(mA/cm2)
		g		(mho/cm2)
		ninf
		ntau 	(ms)
		vhalf 	(mV)
		sf1 	(mV)
		pCa
		cai		(mM)
}

STATE {
        n
}
 
BREAKPOINT {
        SOLVE states METHOD cnexp
        g = gbar*n
		ik = g*(v - ek)
}
 
 
INITIAL {
		rates(v)
		n = ninf
}

DERIVATIVE states {  
        rates(v)
        n' = (ninf-n)/ntau
}
 

PROCEDURE rates(v(mV)) { 
UNITSOFF
		:"n" potassium activation system
		:Data Fit: Scholz et al., 1998
		pCa = log10(cai*1e-3)
		vhalf = (-43.4)*pCa + (-203)
		sf1 =  33.88*exp(-((pCa+5.423)/1.852)^2)
        ninf = 1/(1+exp((vhalf-v)/sf1))
        
		ntau = 5.54522*exp(-v/-42.90548)+0.74926-0.11573*v	: Zhang et al., 2010
}
UNITSON
