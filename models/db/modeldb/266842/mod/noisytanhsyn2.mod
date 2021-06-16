TITLE Continuous synaptic tanh-type transmission mechanism
COMMENT
	This mechanism is the same as tanhsyn, but the synaptic transmission delay
	can be specified via the integration loop, i.e., proc advance(), in the main hoc
	file. Specifically, a vector records the membrane voltage traces of the
	presynaptic neurons, and the membrane voltage at a certain previous time
	point can be assigned to vpre in this mechanism.
ENDCOMMENT

NEURON {
	POINT_PROCESS tanhSyn2
	RANGE tau, e, i, noise, alpha, s, g, rparam, voff, vpre
	NONSPECIFIC_CURRENT i
}

UNITS {
	(uS) = (microsiemens)
	(mS) = (millisiemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER {
	e = 0	(mV)
	alpha = 0.55
	tau = 5.26 (ms) <1e-9,1e9>
	noise = 0
	g = 2.7e-4 (uS)
	rparam = 4
	voff = 0 (mV)
	vpre = -60 (mV)
}

ASSIGNED {
	v (mV)
	i (nA)
}

STATE {
	s
}

INITIAL {
	s=0
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	i = g*s*(v - e)
}

DERIVATIVE state {
	s' = alpha*(1+tanh(vpre-voff/rparam))*(1-s) - s/tau + noise
}