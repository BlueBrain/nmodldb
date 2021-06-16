TITLE Store operated calcium channel for bladder small DRG neuron soma model

: Author: Darshan Mandge (darshanmandge@iitb.ac.in)
: Computational Neurophysiology Lab
: Indian Institute of Technology Bombay, India 

: For details refer: 
: A biophysically detailed computational model of bladder small DRG neuron soma 
: Darshan Mandge and Rohit Manchanda, PLOS Computational Biology (2018)

NEURON {
	SUFFIX soce
	USEION ca READ cai, cao WRITE ica
	USEION caer READ caeri VALENCE 2 	
	RANGE minf, mtau
	RANGE pmax, ica
	THREADSAFE
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
	
	pmax = 1e-9	(cm/s)
	
	: Luik et al., 2008
	kd = 169e-3 (mM)  : dissociation constant
	nh = 4.2  (1)     : hill coefficient
}


STATE {
	m 
}

ASSIGNED {
	ica	 (mA/cm2)
	mtau (ms)
	minf (1)
	cai  (mM)
	cao	 (mM)
	caeri (mM)
}

BREAKPOINT { 
	SOLVE state METHOD cnexp
	ica = pmax*m*ghk(v,cai,cao,2)
}

DERIVATIVE state {
	rates(caeri)
	m'= (minf-m) / mtau
}

INITIAL {
	rates(caeri)
	m = minf
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

PROCEDURE rates(c (mM)) {
	UNITSOFF
	minf = 1 - (c^nh)/(c^nh+kd^nh) :reverse hill equation : Luik et al., 2008
    UNITSON
	mtau = 5000 (ms)  
}

UNITSON 
