TITLE Potassium Density Mechanism

UNITS {
	(mV) = (millivolts)
	(mA) = (milliamp)
	(S) = (siemens)
}

NEURON {
	SUFFIX k
	USEION k WRITE ik
	RANGE gmax, e
}

PARAMETER {
	gmax = 0.000373333 (S/cm2)
	e = -85 (millivolt)
}

ASSIGNED {
	v (mV)
	ek (mV)
	ik (mA/cm2)
	g (S/cm2)
}

STATE { n }

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gmax * n^4
	ik = g * (v - e)
}

INITIAL {
	n = alpha(v)/(alpha(v) + beta(v))
}

DERIVATIVE states {
	n' = alpha(v)*(1-n) - beta(v)*n
}

FUNCTION alpha(Vm (mV)) (/ms) {
	UNITSOFF
	alpha = 0.05 * exp(0.125 * (Vm + 29))
	UNITSON
}

FUNCTION beta(Vm (mV)) (/ms) {
	UNITSOFF
	beta = 0.05 * exp(-0.125 * (Vm + 29))
	UNITSON
}
