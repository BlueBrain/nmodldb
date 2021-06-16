COMMENT
This file, bkfast.mod, implements the fast activating potassium g^B_K(fast) 
current from Quadroni and Knopfel 1994 table 1 that is in the Type B cell model
ENDCOMMENT

NEURON {
	SUFFIX bkfast
	NONSPECIFIC_CURRENT i
	RANGE i, Erev, gbar
	:	GLOBAL taun_min
}

UNITS {
	(S)	=	(siemens)
	(mV)	=	(millivolt)
	(mA)	=	(milliamp)
}

PARAMETER {
	gbar = 37530e-6	(S/cm2) < 0, 1e9 >
	Erev = -82 (mV)
	: taun_min = 0.8 (ms)
}

ASSIGNED {
	i (mA/cm2)
	v (mV)
	g (S/cm2)
	ninf
	tau_n (ms)
}

STATE {	n }

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gbar * n*n
	i = g * (v - Erev)
}

INITIAL {
	: assume that v has been constant for a long time
	n = alphan(v)/(alphan(v) + betan(v))
}

DERIVATIVE states {
	rates(v)
	n' = (ninf - n)/tau_n
}

FUNCTION alphan(Vm (mV)) (/ms) {
	UNITSOFF
	alphan = 5.82 /(1 + exp( -0.125 * (Vm + 3.3)))
	UNITSON
}

FUNCTION betan(Vm (mV)) (/ms) {
	UNITSOFF
	betan =  2.413 / (1 + exp( 0.0675 * (Vm + 46.35)))
	UNITSON
}

FUNCTION taun(Vm (mV)) (/ms) {
	UNITSOFF
	taun = 1.0 / (alphan(Vm) + betan(Vm))
	: if (taun < taun_min) {
	:	taun = taun_min
	: }
	UNITSON
}

PROCEDURE rates(Vm(mV)) {
	tau_n = taun(Vm)
	: ninf = alphan(Vm)/(alphan(Vm) + betan(Vm))
	ninf = alphan(Vm) * tau_n	: change back to above if use taun_min
}
