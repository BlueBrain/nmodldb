:
: USAGE
: objref modyn
: elec[0] modyn = new mdltrdyn(0.5)
: dend[0] new NetCon(&stimon_tcihshift, modyn, 0.5, 0, weight)
: 
: setpointer 
: soma[0].pmodyn_tcihshift
:

NEURON {
POINT_PROCESS mdltrdyn
RANGE tau1, tau2
NONSPECIFIC_CURRENT i
}

PARAMETER {
tau1 = 10 (ms)
tau2 = 20 (ms)
}

ASSIGNED {
v (millivolt)
i (nanoamp)
}

STATE { 
	g (nanoamp)
	a (nanoamp)
	b (nanoamp)
}

INITIAL { 
	g = 0
	a = 0
	b = 0
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	i = g
}

DERIVATIVE state {
	a' = -a/tau1
	b' = -b/tau2
 	g = b - a        
}

NET_RECEIVE(weight (nanoamp)) {
	a = a + weight
	b = b + weight
}
