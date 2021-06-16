TITLE L-type calcium channel with high threshold for activation
: used in somatic and dendritic regions 
: 
: After Borg 
:  Updated by Maria Markaki  12/02/03

NEURON {
	SUFFIX cal
	USEION ca READ cai, eca WRITE ica
 	USEION cal WRITE ical VALENCE 2: to use it in cal-specific pump
        RANGE gcalbar, ica, po, ical
	GLOBAL inf, s_inf, tau_m
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) =	(millimolar)
	FARADAY = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
}


PARAMETER {     
  	ki     = 0.025  (mM)            : middle point of inactivation fct
  :	ki     = 0.01  (mM)            : middle point of inactivation fct
 :   gcalbar = 0.2e-7      (cm/s) : initialized conductance
	gcalbar = 0   (mho/cm2)  : initialized conductance
 	taumin  = 100    (ms)            : minimal value of the time cst
        vhalf = -1 (mV)       :half potential for activation 
	zeta=-4.6
	t0=1.5(ms)
	b = 0.01 	(mM)
:	b = 0.005 	(mM)
}


ASSIGNED {      : parameters needed to solve DE
        v               (mV)
 	celsius         (degC)
	cai             (mM)      : initial internal Ca++ concentration
	ica             (mA/cm2)
	eca             (mV)
	ical             (mA/cm2)
	po
        inf
	s_inf
	tau_m           (ms)
}

STATE {	
	m 
	s 
} 


INITIAL {
	rates(v,cai)
	m = inf    : initial activation parameter value
	s = s_inf
}

FUNCTION h2(cai(mM)) {
	h2 = ki/(ki+cai)
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	po = m*m*h2(cai)
:	po = m*m*h2(cai)
	ica = (gcalbar*po+s*s*8*gcalbar)*(v-eca)
:	ica = (gcalbar*po+s*gcalbar*0.4)*ghk(v, cai, cao)
	ical = ica
}


DERIVATIVE states {
	rates(v,cai)
	m' = (inf-m)/t0
	s' = (s_inf-s)/tau_m
}



FUNCTION alp(v(mV)) {       
UNITSOFF
  alp = exp(1.e-3*zeta*(v-vhalf)*9.648e4/(8.315*(273.16+celsius))) 
UNITSON
}

PROCEDURE rates(v(mV), cai(mM)) {LOCAL a, alpha2
		a = alp(v)
		inf = 1/(1+a)
		alpha2 = (cai/b)^2
		s_inf = alpha2 / (alpha2 + 1)
		tau_m = taumin+ 1(ms)*1(mM)/(cai+b)
}

