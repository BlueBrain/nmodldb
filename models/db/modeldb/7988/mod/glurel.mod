
TITLE transmitter release


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX glurel
	RANGE dur, cmax, T, Twait
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mM) = (milli/liter)
}

PARAMETER {
	dur = 2000 (ms)
	cmax = 1 (mM)
	Twait = 10 (ms)
}

ASSIGNED {
	T (mM)
}


INITIAL {
	T = 0
}

BREAKPOINT {
        if(t < Twait) {
                T = 0
        }
	
        if ( (t >= Twait)&&(t <= (Twait + dur))) {
		T = cmax
        }
	  if (t > (Twait + dur)) {
		T = 0
	  }

}

