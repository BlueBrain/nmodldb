: relative refractory period for a nonlinear IF cell
: implemented as a high conductance
: Usage is to have a one of the NetCon's watching the cell voltage
: send an event with 0 delay to this object

NEURON {
	POINT_PROCESS Refrac_rel
	RANGE gr, e,  g, tau
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
}

PARAMETER {
	tau=1 (ms)   : the decay 
	gr = 1 (umho) : conductance during refractory period
	e = -65 (mV) : refractory channel reversal potential
}

ASSIGNED {
	
	v (mV)
	i (nA)
}

STATE {
	g (umho)
}


INITIAL {
	g = 0 : not in refractory period
}


BREAKPOINT {
	SOLVE state METHOD cnexp
	i = g*(v - e)
}

DERIVATIVE state {
	g' = -g/tau
}


NET_RECEIVE(w) { : w not used. external event is for entering refractory period
	if (flag == 0) { : external event
		g = g+gr   
	}
}

