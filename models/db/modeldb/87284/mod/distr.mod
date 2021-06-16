TITLE ...just to store peak membrane voltage
: M.Migliore June 2001
: T Morse February 2010 added times of occurence

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {
	v (mV)
}


NEURON {
	SUFFIX ds
        RANGE vmax, tmax
}

ASSIGNED {
	vmax
	tmax
}

INITIAL {
	vmax=v
}


BREAKPOINT {
	if (v>vmax) {vmax=v tmax=t}
}
