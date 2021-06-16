COMMENT
This file, ahp.mod, implements the after hyperpolarization (gAHP)  
current from Quadroni and Knopfel 1994 table 1
Note: this channel can be verified by testing if gives behavior specified in Quadroni
thesis, i.e. that channel is inactivated at cai concentration of 50 nM (5e-5 mM) and 
essentially activated at 500 nM (5e-4 mM) (at ten fold higher concentration )
ENDCOMMENT

NEURON {
	SUFFIX ahp
	USEION ca READ cai
	NONSPECIFIC_CURRENT i
	RANGE i, Erev, gbar, q, tau_q, qinf, betaq_const
}

UNITS {
	(S)	=	(siemens)
	(mV)	=	(millivolt)
	(mA)	=	(milliamp)
}

PARAMETER {
	gbar = 2167e-6	(S/cm2) < 0, 1e9 >
	Erev = -82 (mV)
	cai (mM) : starts at 0.050 uM = 5e-8 (M) = 5e-5 mM
	betaq_const =  0.074	: Q+K 94
}

ASSIGNED {
	i (mA/cm2)
	v (mV)
	g (S/cm2)
	qinf
	tau_q (ms)
}

STATE {	q }

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gbar * q*q
	i = g * (v - Erev)
}

INITIAL {
	: assume that v has been constant for a long time
	q = alphaq(v)/(alphaq(v) + betaq(v))
}

DERIVATIVE states {
	rates(v)
	q' = (qinf - q)/tau_q
}

FUNCTION alphaq(Vm (mV)) (/ms) {
	UNITSOFF
	alphaq = 3.5e9 * (cai)^3
	UNITSON
}

FUNCTION betaq(Vm (mV)) (/ms) {
	UNITSOFF
	betaq = betaq_const
	UNITSON
}

FUNCTION tauq(Vm (mV)) (/ms) {
	UNITSOFF
	tauq = 1.0 / (alphaq(Vm) + betaq(Vm))
	UNITSON
}

PROCEDURE rates(Vm(mV)) {
	tau_q = tauq(Vm)
	qinf = alphaq(Vm) * tau_q      : change back to a/(a+b) if use tauq_min
}
