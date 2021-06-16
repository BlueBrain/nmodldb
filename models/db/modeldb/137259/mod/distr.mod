TITLE ...just to store peak membrane voltage
: M.Migliore June 2001
: J Baker added del, dur, vmin Jan 2008

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	v (mV)
	del=0 (ms)
	dur=1e9 (ms)
}


NEURON {
	SUFFIX ds
	RANGE vmax, vmin
}

ASSIGNED {
	vmax (mV)
	vmin (mV)
}

INITIAL {
	vmax=-1e9
	vmin=1e9
}


BREAKPOINT {
	if (t>=del && t<=del+dur) {
		if (v>vmax) {vmax=v}
		if (v<vmin) {vmin=v}
	}
}
