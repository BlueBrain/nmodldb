COMMENT
This file, a.mod, implements the transient potassium current gA from 
Quadroni and Knopfel 1994 table 1
ENDCOMMENT

NEURON {
	SUFFIX a
	NONSPECIFIC_CURRENT i
	RANGE i, Erev, gbar, a, b
	GLOBAL taua_min, taub_min
}

UNITS {
	(S)	=	(siemens)
	(mV)	=	(millivolt)
	(mA)	=	(milliamp)
}

PARAMETER {
	gbar = 1829e-6	(S/cm2) < 0, 1e9 >
	Erev = -82 (mV)
	taua_min = 1.0 (ms)
	taub_min = 24.0 (ms)
}

ASSIGNED {
	i (mA/cm2)
	v (mV)
	g (S/cm2)
	ainf
	binf
	tau_a (ms)
	tau_b (ms)
}

STATE {	a b }

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gbar * a^3 * b
	i = g * (v - Erev)
}

INITIAL {
	: assume that v has been constant for a long time
	a = alphaa(v)/(alphaa(v) + betaa(v))
	b = alphab(v)/(alphab(v) + betab(v))
}

DERIVATIVE states {
	rates(v)
	a' = (ainf - a)/tau_a
	b' = (binf - b)/tau_b
}

FUNCTION alphaa(Vm (mV)) (/ms) {
	UNITSOFF
	alphaa = 0.2 * exp( 0.14 * (Vm + 65.0))
	UNITSON
}

FUNCTION betaa(Vm (mV)) (/ms) {
	UNITSOFF
	betaa =  0.2 * exp( -0.035 * (Vm + 65.0))
	UNITSON
}

FUNCTION taua(Vm (mV)) (/ms) {
	UNITSOFF
	taua = 1.0 / (alphaa(Vm) + betaa(Vm))
	if (taua < taua_min) {
		taua = taua_min
	}
	UNITSON
}

FUNCTION alphab(Vm (mV)) (/ms) {
	UNITSOFF
	alphab = 0.01 * exp( -0.11 * (Vm + 71.0))
	UNITSON
}

FUNCTION betab(Vm (mV)) (/ms) {
	UNITSOFF
	betab =  0.01 * exp( 0.164 * (Vm + 71.0))
	UNITSON
}

FUNCTION taub(Vm (mV)) (/ms) {
	UNITSOFF
	taub = 1.0 / (alphab(Vm) + betab(Vm))
	if (taub < taub_min) {
		taub = taub_min
	}
	UNITSON
}

PROCEDURE rates(Vm(mV)) {
	tau_a = taua(Vm)
	ainf = alphaa(Vm)/(alphaa(Vm) + betaa(Vm))
	tau_b = taub(Vm)
	binf = alphab(Vm)/(alphab(Vm) + betab(Vm))
}
