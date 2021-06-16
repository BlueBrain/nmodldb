COMMENT
** This mod file is based on the alpha synapse in NEURON, except for a ratio factor con in the NET_RECEIVE block.**
ENDCOMMENT

NEURON {
	POINT_PROCESS Alpha
	RANGE tau1, tau2, e, i,con
	NONSPECIFIC_CURRENT i

	RANGE g
	GLOBAL total
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
}

PARAMETER {
	tau1=.1 (ms) <1e-9,1e9>
	tau2 = 10 (ms) <1e-9,1e9>
	e=0	(mV)
        con=10
}

ASSIGNED {
	v (mV)
	i (nA)
	g (umho)
	factor
	total (umho)
}

STATE {
	A (umho)
	B (umho)
}

INITIAL {
	LOCAL tp
	total = 0
	if (tau1/tau2 > .9999) {
		tau1 = .9999*tau2
	}
	A = 0
	B = 0
	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
	factor = -exp(-tp/tau1) + exp(-tp/tau2)
	factor = 1/factor
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	g = B - A
	i = g*(v - e)
}

DERIVATIVE state {
	A' = -A/tau1
	B' = -B/tau2
}

NET_RECEIVE(weight (umho)) {
	state_discontinuity(A, A + weight*factor*con)
	state_discontinuity(B, B + weight*factor*con)
	total = total+weight*con
}
