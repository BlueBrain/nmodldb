TITLE L-type HVA calcium current for bladder small DRG neuron soma model

: Author: Darshan Mandge (darshanmandge@iitb.ac.in)
: Computational Neurophysiology Lab
: Indian Institute of Technology Bombay, India 

: For details refer: 
: A biophysically detailed computational model of bladder small DRG neuron soma 
: Darshan Mandge and Rohit Manchanda, PLOS Computational Biology (2018)

NEURON {
	SUFFIX CaL
	USEION ca READ cai,cao WRITE ica
	RANGE minf, mtau, hinf, htau, ica
	RANGE pmax, hca
}

UNITS {
	(mA)	= (milliamp)
	(mV)	= (millivolt)
	(mM)	= (milli/liter)
    FARADAY = 96489 (coul)
     R      = 8.314 (volt-coul/degC)
}

PARAMETER {
    v		(mV)
    celsius	(degC)
    cao		(mM)
	
    : Calcium dependent inactivation
    cai0 = 136e-6	(mM)
    Kd	 = 1e-3 	(mM)
	pmax = 2.7e-5	(cm/s)	
	hca 			(1)
}

STATE {
	m h
}

ASSIGNED {
	ica			(mA/cm2)
	cai			(mM)
    mtau		(ms)
	minf
	hinf
	htau		(ms)
}

BREAKPOINT { 
	SOLVE state METHOD cnexp

    hca = 1/(1+(cai/Kd)^4) : Tong et al.,2011 
	ica = pmax*m*h*hca*ghk(v,cai,cao,2)
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
    hca = 1/(1+(cai/Kd)^4)
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

PROCEDURE rates(v(mV)) {
	UNITSOFF
	mtau = 2.10686 + (77.425/(16.02197*sqrt(3.14/2)))*exp(-2*((v+10)/16.022(mV))^2) :Data Fit: Fox et al., 1987
	htau = 825.80 + (31780.18/(39.75*sqrt(3.14/2)))*exp(-2*((v-0)/39.74938(mV))^2) :Data Fit: Fox et al., 1987
	

	minf = 1 / (1+exp((8.45687-v)/(4.26269)))  :Data Fit: Fox et al., 1987
	hinf = 1 / (1+exp((v+42.52246)/(7.47793))) :Data Fit: Fox et al., 1987
	UNITSON
}

