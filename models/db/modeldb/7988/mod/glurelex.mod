TITLE transmitter release

: exponential increase and decrease of glutamate concentration, instead
: of square pulse as in glurel.mod.
: Modified from glurel.mod.


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX glurelex
	RANGE dur, cmax, T, Twait, tau1, tau2, tphase2, cmaxphase2
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mM) = (milli/liter)
}

PARAMETER {
	dur = 1 (ms)
	cmax = 2.7 (mM)
	Twait = 10 (ms)
	cmaxphase2 = 0.41 (mM)
	tau1 = 0.1 (ms)
	tau2 = 2.1 (ms)
}

ASSIGNED {
	T (mM)
	tphase2 (ms)
}


INITIAL {
	T = 0
}


BREAKPOINT {
        if(t < Twait) {
                T = 0
        }
	
	  if(t == Twait) {
		 T = cmax + cmaxphase2
	  }

	  if ( t >= Twait ) {
		if ( T == 0) {
			T = cmax + cmaxphase2
		}
	  	T = cmax*exp(-(t - Twait)/tau1) + cmaxphase2*exp(-(t - Twait)/tau2)	
	}	

}

