TITLE Slow Ca-dependent potassium current
:
:   Ca++ dependent K+ current IC responsible for slow AHP

NEURON {
	SUFFIX kca
	USEION k READ ek WRITE ik
	USEION cal READ cali
	RANGE  gbar, po, ik
	GLOBAL m_inf, tau_m
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
}

ASSIGNED {       : parameters needed to solve DE
	v               (mV)
	celsius         (degC)
	ek              (mV)
	cali             (mM)           : initial [Ca]i
	ik              (mA/cm2)
	po
	m_inf
	tau_m           (ms)
}

PARAMETER {
	gbar    = 0.01   (mho/cm2)
 	taumin  = 100    (ms)            : minimal value of the time cst
:	b = 0.1		(/ms)
	b = 0.005 	(mM)
:	a0 = 250 (/ms-mM)
}


STATE {
	m
}

BREAKPOINT { 
	SOLVE states METHOD cnexp
:	po = m
	po = m*m
	ik = gbar*po*(v + 85)    : potassium current induced by this channel
:	ik = gbar*po*(v - ek)    : potassium current induced by this channel
}

DERIVATIVE states {
	rates(cali)
	m' = (m_inf - m) / tau_m
}


INITIAL {
	rates(cali)
	m = m_inf
}

:FUNCTION alp(cai (mM))(/ms){
:	alp = a0*cai
:	}

PROCEDURE rates(cali(mM)) {  LOCAL a 
:	a = alp(cali)
	a = cali/b
 	m_inf = a/(a+1)
	tau_m = taumin+ 1(ms)*1(mM)/(cali+b)
}
