COMMENT
Implements the "spike" mechanism used by Gruber et al. 2003.
Watch this mechanism's spike variable with a NetCon that has 
threshold 0.5, and record that NetCon's spike event times.

NOTES:
  Needs fixed time step integration, 
    won't work with adaptive integration.
  Needs NEURON v. 5.5 or later.
ENDCOMMENT

NEURON {
	POINT_PROCESS GSpike
	RANGE vt, vh, vc, spike, minisi
}

UNITS {
	(mV) = (millivolt)
}

PARAMETER {
	vt = -58 (mV)
	vh = -55 (mV)
	vc = 2.5 (mV)
	minisi = 20 (ms)
	spikewidth = 1 (ms)
}

ASSIGNED {
	v (mV)
	tp (ms)
	spike (1)
}

INITIAL {
	tp = -1e9 (ms)
	spike = 0
}

AFTER SOLVE {
	if (at_time(tp+spikewidth)) {
		spike = 0
	}
	if (v>vt) {
		if (t>tp+minisi*(1+exp(-(v-vh)/vc))) {
			spike = 1
			tp = t
		}
	}
}
