TITLE McCormick and Huguenard Na and K leak channels

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX tcleakdepol
	NONSPECIFIC_CURRENT inal
	USEION k READ ek WRITE ik
	RANGE gnal, gkl, ena, minc , tau_m, stimon, nsm, iksub
	POINTER vext, pmodyn		: vext turns into a dummy variable to determine whether a stim pulse is ON
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gnal	= 0.00000086	(mho/cm2)
	gkl  	= 0.03	(mho/cm2)
	m_inf = 1 
	mmin = 0
	tau_m = 2000 (ms)
	minc = 0.0025 (/ms)
	stimon = 0
	nsm = 1

:	ek	= -95	(mV)
	ena	= 45	(mV)

	celsius		(degC)
	dt              (ms)
	v               (mV)
}

STATE {
	m
}

ASSIGNED {
	ek (mV)
	inal	(mA/cm2)
	iksub	(mA/cm2)
	ik	(mA/cm2)
	vext	: vext turns into a dummy variable to determine whether a stim pulse is ON
	pmodyn
}


BREAKPOINT {
	SOLVE states METHOD euler
	
	inal = gnal * (v - ena)
	
	if (m < 0) {
		m = 0
	}
	iksub = gkl * m * (v - ek)
	ik = iksub
}

DERIVATIVE states {   : exact Hodgkin-Huxley equations
	eval_stimon(vext)
	m' = (m_inf - m) / tau_m - nsm*pmodyn*minc
}

INITIAL {
	m=m_inf
}

PROCEDURE eval_stimon(vext) {  : vext turns into a dummy variable to determine whether a stim pulse is ON
	if(vext==0) {
		stimon=0
	}else{
		stimon=1
	}
}