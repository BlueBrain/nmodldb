TITLE Continuous and instantaneous tanh-type synaptic transmission mechanism
COMMENT
	The presynaptic membrane potential can be specified via setpointer
	in the hoc file.
ENDCOMMENT

NEURON {
	POINT_PROCESS tanhSyn
	RANGE tau, e, i, noise, alpha, s, g, rparam, voff
	NONSPECIFIC_CURRENT i
	POINTER vpre
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
}

ASSIGNED {
	v (mV)
	vpre (mV)
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