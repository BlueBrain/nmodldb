TITLE Cardiac time independent inward rectifier IK1 current
: from BEELER & REUTER, J.Physiol, 1977

NEURON {
	SUFFIX IK1
	USEION k WRITE ik
	RANGE gK1, ik
	GLOBAL dummy : prevent vectorization for use with CVODE
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mM) = (milli/liter)
	(S) = (siemens)
}

PARAMETER { 
	gK1=0.00035 (S/cm2) <0,1e9>
}

ASSIGNED {
	v (mV)
	ik (mA/cm2)
	dummy
}

BREAKPOINT {
LOCAL d, n, r
	d = 4 * (exp(0.04*(v + 85))-1)
	n = (exp(0.08*(v + 53)) + exp(0.04*(v + 53)))
	r = 0.2*(v + 23)/(1-exp(-0.04*(v + 23)))
	ik = gK1*(d/n + r)
}