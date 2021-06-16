: step pulse with ramp increase and ramp decrease.



TITLE transmitter release


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX relramp
	RANGE dur, cmax, T, Twait, rampdur
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mM) = (milli/liter)
}

PARAMETER {
	dur = 20 (ms)
	cmax = 1 (mM)
	Twait = 10 (ms)
	rampdur = 4 (ms)
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
	
        if ( (t >= Twait)&&(t <= (Twait + rampdur))) {
		T = cmax*(t - Twait)/rampdur
        }
        if ((t > ( Twait + rampdur)) && (t <= (Twait + dur))) {

		T = cmax
        }
	  if ((t > (Twait + dur)) &&	(t <= (Twait + dur + rampdur))) {
		T = cmax - cmax*(t - Twait - dur)/rampdur
	  }
	  if (t > (Twait + dur + rampdur)) {
		T = 0
	  }

}

