TITLE R-type calcium current for bladder small DRG neuron soma model

: Author: Darshan Mandge (darshanmandge@iitb.ac.in)
: Computational Neurophysiology Lab
: Indian Institute of Technology Bombay, India 

: For details refer: 
: A biophysically detailed computational model of bladder small DRG neuron soma 
: Darshan Mandge and Rohit Manchanda, PLOS Computational Biology (2018)

NEURON {
	SUFFIX CaR
	USEION ca READ cai,cao WRITE ica
	RANGE minf, mtau, hinf, htau, h2tau, ica, pmax
}

UNITS {
	(mA)	= (milliamp)
	(mV)	= (millivolt)
	(mM)	= (milli/liter)
	FARADAY = 96489 (coul)
	R       = 8.314 (volt-coul/degC)
}

PARAMETER {
	v		(mV)
	celsius	(degC)
	cai		(mM)
	cao		(mM)
	pmax =  1e-8	(cm/s)
}

STATE {
	m 	:	activation
	h 	:	slow inactivation
	h2	: 	fast inactivation
}

ASSIGNED {
	ica		(mA/cm2)
	mtau		(ms)
	minf
	hinf
	htau		(ms)
    h2tau		(ms)
}

BREAKPOINT { 
	SOLVE state METHOD cnexp
	ica = pmax*m*(0.4*h+0.6*h2)*ghk(v,cai,cao,2)
}

DERIVATIVE state {
	rates(v)
	m'= (minf-m) / mtau
	h'= (hinf-h) / htau
    h2'= (hinf-h2) / h2tau
}

INITIAL {
	rates(v)
	m = minf
	h = hinf
    h2= hinf
}

: For GHK Current equation (adapted from Stepheen and Manchanda, 2009)
FUNCTION ghk( v(mV), ci(mM), co(mM), z)  (millicoul/cm3) { LOCAL e, w
        w = v * (.001) * z*FARADAY / (R*(celsius+273.16))
        if (fabs(w)>1e-4) 
          { e = w / (exp(w)-1) }
        else : denominator is small -> Taylor series
          { e = 1-w/2 }
        ghk = - (.001) * z*FARADAY * (co-ci*exp(w)) * e
}

UNITSOFF

FUNCTION_TABLE tabmtau(v(mV)) (ms) 
FUNCTION_TABLE tabhtau(v(mV)) (ms)
FUNCTION_TABLE tabhtau2(v(mV)) (ms)


PROCEDURE rates(v(mV)) { 	
	minf = 1 / (1+exp((-5-v)/5)) 	: Li et al., 2007
	mtau = tabmtau(v)				: Hilaire et al., 1997
	
	hinf = 1 / (1+exp((v+51)/12)) 	: Hilaire et al., 1997
	    
    htau  = tabhtau(v)	: Hilaire et al., 1997
	h2tau = tabhtau2(v) : Hilaire et al., 1997
}

UNITSON 
