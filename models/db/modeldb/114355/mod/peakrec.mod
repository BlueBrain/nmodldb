COMMENT
record the peak of voltage in a compartment
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

DEFINE SIZE 100

NEURON {
	POINT_PROCESS PeakRec
	RANGE time,peak
}

UNITS {
	(mV) = (millivolt)
}

ASSIGNED {
	peak
	time
	v_init
}

INITIAL {
	v_init = v
	peak = 0
	time = t
}

BREAKPOINT {
	SOLVE check
}

PROCEDURE check() {
VERBATIM
	if (v > peak+v_init ) {
		peak = v-v_init;
		time = t;
	}
	return 0;
ENDVERBATIM
}
