COMMENT
This file, akfast.mod, implements the fast activating potassium g^A_K(fast) 
current from Quadroni and Knopfel 1994 table 1 that is in the Type A cell model
ENDCOMMENT

NEURON {
	SUFFIX akfast
	NONSPECIFIC_CURRENT i
	RANGE i, Erev, gbar, n, tau_n, ninf
	GLOBAL taun_min
}

UNITS {
	(S)	=	(siemens)
	(mV)	=	(millivolt)
	(mA)	=	(milliamp)
}

PARAMETER {
	gbar = 24073e-6	(S/cm2) < 0, 1e9 >
	Erev = -82 (mV)
	taun_min = 0.8 (ms)
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
	g = gbar * n^3
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
	alphan = 0.16 * exp( 0.185 * (Vm + 42.3))
	UNITSON
}

FUNCTION betan(Vm (mV)) (/ms) {
	UNITSOFF
	betan =  0.16 * exp( -0.033 * (Vm + 42.3))
	UNITSON
}

FUNCTION taun(Vm (mV)) (/ms) {
	UNITSOFF
	taun = 1.0 / (alphan(Vm) + betan(Vm))
	if (taun < taun_min) {
		taun = taun_min
	}
	UNITSON
}

PROCEDURE rates(Vm(mV)) {
	tau_n = taun(Vm)
	ninf = alphan(Vm)/(alphan(Vm) + betan(Vm))
}
