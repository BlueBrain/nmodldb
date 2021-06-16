: relative refractory period for a nonlinear IF cell
: implemented as a high conductance
: Usage is to have a one of the NetCon's watching the cell voltage
: send an event with 0 delay to this object

NEURON {
	POINT_PROCESS AHP
	USEION k READ ek 
	RANGE gr, e,  g, tau
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(uS)   = (microsiemens)
}

PARAMETER {
	tau=5 (ms)   : the decay 
	gr = 300e-6 (uS) : conductance during refractory period
	e = -90 (mV) : refractory channel reversal potential
}

ASSIGNED {
	
	v (mV)
	i (nA)
	ek (mV)
}

STATE {
	g (uS)
}


INITIAL {
	g = 0 : not in refractory period
}


BREAKPOINT {
	SOLVE state METHOD cnexp
	i = g*(v - ek)
}

DERIVATIVE state {
	g' = -g/tau
}


NET_RECEIVE(weight (uS)) { : w not used. external event is for entering refractory period
	if (flag == 0) { : external event
		g = g+weight   
	}
}