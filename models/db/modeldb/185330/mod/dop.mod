TITLE simplified Dopamine dynamics

: Modified from the original AMPA.mod, M.Migliore Nov 2007


NEURON {
	POINT_PROCESS dopnet
	RANGE vmax, km
	NONSPECIFIC_CURRENT i
}

UNITS {
}

PARAMETER {

	vmax	= 5e-3		
	km	= 0.2		
}


ASSIGNED {
	i
}

STATE {dop}

INITIAL {
	dop=0
}

BREAKPOINT {
	SOLVE release METHOD cnexp
	i=0
}

DERIVATIVE release {
	dop' = -vmax/((km/dop)+1)
}

NET_RECEIVE(weight) {
	state_discontinuity(dop, dop + weight)
}

