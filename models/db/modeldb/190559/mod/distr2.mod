TITLE ...just to store peak membrane voltage
: M.Migliore June 2001

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {
	v (mV)
}


NEURON {
	SUFFIX ds2
        RANGE vmax1, vmax2
}

ASSIGNED {
	vmax
	vmax1
	vmax2
}

INITIAL {
	vmax1=v
	vmax2=v
}


BREAKPOINT {
	if (t<1000 && v>vmax1) {vmax1=v} else {if (t>34000 && v>vmax2) {vmax2=v}}
}
