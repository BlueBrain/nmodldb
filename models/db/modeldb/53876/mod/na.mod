COMMENT
This file, na.mod, implements the gna current from 
Quadroni and Knopfel 1994 table 1
ENDCOMMENT

NEURON {
	SUFFIX na
	NONSPECIFIC_CURRENT i
	RANGE i, Erev, gbar
	GLOBAL taum_min, tauh_min
}

UNITS {
	(S)	=	(siemens)
	(mV)	=	(millivolt)
	(mA)	=	(milliamp)
}

PARAMETER {
	gbar = 20172e-6	(S/cm2) < 0, 1e9 >
	Erev = 50 (mV)
	taum_min = 0.05 (ms)
	tauh_min = 0.3 (ms)
}

ASSIGNED {
	i (mA/cm2)
	v (mV)
	g (S/cm2)
	minf
	hinf
	tau_h (ms)
	tau_m (ms)
}

STATE {	m h }

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gbar * m^3 * h^2
	i = g * (v - Erev)
}

INITIAL {
	: assume that v has been constant for a long time
	m = alpham(v)/(alpham(v) + betam(v))
	h = alphah(v)/(alphah(v) + betah(v))
}
DERIVATIVE states {
	rates(v)
	m' = (minf - m)/tau_m
	h' = (hinf - h)/tau_h
}

FUNCTION alpham(Vm (mV)) (/ms) {
	UNITSOFF
	alpham = 5.0 * exp( 0.1 * (Vm + 39.5))
	UNITSON
}

FUNCTION betam(Vm (mV)) (/ms) {
	UNITSOFF
	betam =  5.9 * exp( -0.044 * (Vm + 39.5))
	UNITSON
}

FUNCTION alphah(Vm (mV)) (/ms) {
	UNITSOFF
	alphah = 0.567 * exp( -0.024 * (Vm + 35.0))
	UNITSON
}

FUNCTION betah(Vm (mV)) (/ms) {
	UNITSOFF
	betah = 0.567 * exp(0.275 * (Vm + 35.0))
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

FUNCTION tauh(Vm (mV)) (/ms) {
	UNITSOFF
	tauh = 1.0 / (alphah(Vm) + betah(Vm))
	if (tauh < tauh_min) {
		tauh = tauh_min
	}
	UNITSON
}

PROCEDURE rates(Vm(mV)) {
	tau_h = tauh(Vm)
	tau_m = taum(Vm)
	minf = alpham(Vm)/(alpham(Vm) + betam(Vm))
	hinf = alphah(Vm)/(alphah(Vm) + betah(Vm))
}
