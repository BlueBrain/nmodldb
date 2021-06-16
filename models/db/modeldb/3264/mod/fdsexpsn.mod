COMMENT
Implementation of the model of short-term facilitation and depression described in
  Varela, J.A., Sen, K., Gibson, J., Fost, J., Abbott, L.R., and Nelson, S.B.
  A quantitative description of short-term plasticity at excitatory synapses 
  in layer 2/3 of rat primary visual cortex
  Journal of Neuroscience 17:7926-7940, 1997
This is a modification of ExpSyn that can receive multiple streams of 
synaptic input via NetCon objects.  Each stream keeps track of its own 
weight and activation history.

The printf() statements are for testing purposes only.
ENDCOMMENT


NEURON {
	POINT_PROCESS FDSExpSyn
	RANGE tau, e, i, f, tau_F, d1, tau_D1, d2, tau_D2
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
}

PARAMETER {
	tau = 0.1 (ms) < 1e-9, 1e9 >
	e = 0	(mV)
	: these values are from Fig.3 in Varela et al. 1997
	: the (1) is needed for the range limits to be effective
        f = 0.917 (1) < 0, 1e9 >    : facilitation
        tau_F = 94 (ms) < 1e-9, 1e9 >
        d1 = 0.416 (1) < 0, 1 >     : fast depression
        tau_D1 = 380 (ms) < 1e-9, 1e9 >
        d2 = 0.975 (1) < 0, 1 >     : slow depression
        tau_D2 = 9200 (ms) < 1e-9, 1e9 >
}

ASSIGNED {
	v (mV)
	i (nA)
}

STATE {
	g (umho)
}

INITIAL {
	g=0
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	i = g*(v - e)
}

DERIVATIVE state {
	g' = -g/tau
}

NET_RECEIVE(weight (umho), F, D1, D2, tsyn (ms)) {
INITIAL {
: these are in NET_RECEIVE to be per-stream
	F = 1
	D1 = 1
	D2 = 1
	tsyn = t
: this header will appear once per stream
: printf("t\t t-tsyn\t F\t D1\t D2\t amp\t newF\t newD1\t newD2\n")
}

	F = 1 + (F-1)*exp(-(t - tsyn)/tau_F)
	D1 = 1 - (1-D1)*exp(-(t - tsyn)/tau_D1)
	D2 = 1 - (1-D2)*exp(-(t - tsyn)/tau_D2)
: printf("%g\t%g\t%g\t%g\t%g\t%g", t, t-tsyn, F, D1, D2, weight*F*D1*D2)
	tsyn = t

	state_discontinuity(g, g + weight*F*D1*D2)

	F = F + f
	D1 = D1 * d1
	D2 = D2 * d2
: printf("\t%g\t%g\t%g\n", F, D1, D2)
}
