COMMENT
This file, nmda.mod, implements the NMDA-type glutamate receptors
current from Quadroni and Knopfel 1994 equations (9)
This mechanism is applied uniformly to all dendrites in the above paper.
ENDCOMMENT

NEURON {
	SUFFIX nmda
	NONSPECIFIC_CURRENT i
	RANGE i, Erev, gbar
}

UNITS {
	(S)	=	(siemens)
	(mV)	=	(millivolt)
	(mA)	=	(milliamp)
}

PARAMETER {
	gbar = 2570e-6	(S/cm2) < 0, 1e9 >
	Erev = 0 (mV)
	mg_o = 2 (mM)	: Magnesium concentration outside cell
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
	g = gbar * n
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
	alphan = 3.0 * exp(0.035 * (Vm + 10.0))
	UNITSON
}

FUNCTION betan(Vm (mV)) (/ms) {
	UNITSOFF
	betan = mg_o * exp(-0.035 * (Vm + 10.0))
	UNITSON
}

FUNCTION taun(Vm (mV)) (/ms) {
	UNITSOFF
	taun = 1.0 / (alphan(Vm) + betan(Vm))
	UNITSON
}

PROCEDURE rates(Vm(mV)) {
	tau_n = taun(Vm)
	ninf = alphan(Vm) * tau_n      : change back to a/(a+b) if use taun_min
}
