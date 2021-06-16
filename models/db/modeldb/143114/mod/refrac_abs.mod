: absolute refractory period for a nonlinear IF cell
: implemented as a high conductance
: Usage is to have a one of the NetCon's watching the cell voltage
: send an event with 0 delay to this object

NEURON {
	POINT_PROCESS Refrac_abs
	RANGE gr, e, tr, g
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
}

PARAMETER {
	tr = 2 (ms) : refractory period
	gr = 1 (umho) : conductance during refractory period
	e = -65 (mV) : refractory channel reversal potential
}

ASSIGNED {
	g (umho)
	v (mV)
	i (nA)
}

INITIAL {
	g = 0 : not in refractory period
}


BREAKPOINT {
	i = g*(v - e)
}

NET_RECEIVE(w) { : w not used. external event is for entering refractory period
	if (flag == 0) { : external event
		g = gr
		net_send(tr, 1)
	}else{ : refractory period over
		g = 0
	}
}

