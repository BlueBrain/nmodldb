COMMENT
This file, nap.mod, implements the persistent sodium (Gnap) current from 
Quadroni and Knopfel 1994 table 1
ENDCOMMENT

NEURON {
	SUFFIX nap
	NONSPECIFIC_CURRENT i
	RANGE i, Erev, gbar
	GLOBAL taum_min
}

UNITS {
	(S)	=	(siemens)
	(mV)	=	(millivolt)
	(mA)	=	(milliamp)
}

PARAMETER {
	gbar = 33e-6	(S/cm2) < 0, 1e9 >
	Erev = 50 (mV)
	taum_min = 1.0 (ms)
}

ASSIGNED {
	i (mA/cm2)
	v (mV)
	g (S/cm2)
	minf
	tau_m (ms)
}

STATE {	m }

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gbar * m
	i = g * (v - Erev)
}

INITIAL {
	: assume that v has been constant for a long time
	m = alpham(v)/(alpham(v) + betam(v))
}

DERIVATIVE states {
	rates(v)
	m' = (minf - m)/tau_m
}

FUNCTION alpham(Vm (mV)) (/ms) {
	UNITSOFF
	alpham = 0.12 * exp( 0.12 * (Vm + 56))
	UNITSON
}

FUNCTION betam(Vm (mV)) (/ms) {
	UNITSOFF
	betam =  0.12 * exp( -0.03 * (Vm + 56))
	UNITSON
}

FUNCTION taum(Vm (mV)) (/ms) {
	UNITSOFF
	taum = 1.0 / (alpham(Vm) + betam(Vm))
	if (taum < taum_min) {
		taum = taum_min
	}
	UNITSON
}

PROCEDURE rates(Vm(mV)) {
	tau_m = taum(Vm)
	minf = alpham(Vm)/(alpham(Vm) + betam(Vm))
}
