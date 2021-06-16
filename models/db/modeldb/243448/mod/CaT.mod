TITLE T-type calcium current for bladder small DRG neuron soma model

: Author: Darshan Mandge (darshanmandge@iitb.ac.in)
: Computational Neurophysiology Lab
: Indian Institute of Technology Bombay, India 

: For details refer: 
: A biophysically detailed computational model of bladder small DRG neuron soma 
: Darshan Mandge and Rohit Manchanda, PLOS Computational Biology (2018)

NEURON {
	SUFFIX CaT
	USEION ca READ cai, cao WRITE ica VALENCE 2
    RANGE minf, mtau, hinf, htau, ica, pmax
}

UNITS {
	(mA)	= (milliamp)
	(mV)	= (millivolt)
	(mM)	= (milli/liter)
	FARADAY = 96489 (coul)
	R       = 8.314 (volt-coul/degC)
	(um)	= (micron)
}

PARAMETER {
	v		(mV)
	celsius	(degC)
	cai		(mM)
	cao		(mM)
	pmax = 1e-8	(cm/s)
}

STATE {
	m h
}

ASSIGNED {
	ica		(mA/cm2)
	mtau	(ms)
	minf
	hinf
	htau	(ms)
	area	(um2)
}

BREAKPOINT { 
	SOLVE state METHOD cnexp
	ica = pmax *m*h*ghk(v,cai,cao,2)
}

DERIVATIVE state {
	rates(v)
	m'= (minf-m) / mtau
	h'= (hinf-h) / htau
}

INITIAL {
	rates(v)
	m = minf
	h = hinf
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
FUNCTION_TABLE tabmtau(v(mV)) (ms)
FUNCTION_TABLE tabhtau(v(mV)) (ms)
UNITSOFF

PROCEDURE rates(v(mV)) { 
	mtau = tabmtau(v) : Adapted from Steephen and Manchanda, 2009
	htau = tabhtau(v) : Data from Fox et al.,1987

	minf = 1 / (1+exp((-55.29-v)/6.38)) : Data Fit: Fox et al.,1987
	hinf = 1 / (1+exp((v+76.59)/4.46))  : Data Fit: Fox et al.,1987
}
