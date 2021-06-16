TITLE N-type calcium current for bladder small DRG neuron soma model

: Author: Darshan Mandge (darshanmandge@iitb.ac.in)
: Computational Neurophysiology Lab
: Indian Institute of Technology Bombay, India 

: For details refer: 
: A biophysically detailed computational model of bladder small DRG neuron soma 
: Darshan Mandge and Rohit Manchanda, PLOS Computational Biology (2018)

NEURON {
	SUFFIX CaN
	USEION ca READ cai,cao WRITE ica
	RANGE minf, mtau, hinf, htau, ica
	RANGE pmax, a, hca
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
	cao		(mM)
	pmax 	= 2.8e-5 (cm/s)
	a		= 0.7326 (1)
	
    cai0 = 136e-6 (mM) 
    Kd = 1e-3     (mM)
	hca 		  (1)
}

STATE {
	m h
}

ASSIGNED {
	ica		(mA/cm2)
    cai		(mM)
	mtau	(ms)
	minf
	hinf
	htau	(ms)
}

BREAKPOINT { 
	SOLVE state METHOD cnexp
	hca = 1/(1+(cai/Kd)^4) : Tong et al, 2011. Calcium dependent inactivation
	ica = pmax*m*(a*h+(1-a))*hca*ghk(v,cai,cao,2)
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

UNITSOFF
FUNCTION_TABLE tabhtau(v(mV)) (ms)

 PROCEDURE rates(v(mV)) {
	mtau = 0.80448 + (101.05111/(14.9947*sqrt(3.14/2)))*exp(-2*((v+20)/14.9947)^2)
 
	htau = tabhtau(v)  : Aosaki et al., 1989

	minf = 1 / (1+exp((-6.5-v)/6.5)) : Fox et al., 1987
	hinf = 1 / (1+exp((v+70)/12.5))  : Fox et al., 1987
}

UNITSON 