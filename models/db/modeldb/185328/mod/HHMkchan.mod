TITLE HHMk channel
: Hodgkin Huxley Moore k channel Mod


NEURON {
	SUFFIX HHkM
	USEION k READ ek WRITE ik
	RANGE gkbar, ik
	GLOBAL inf
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
        (S) = (siemens)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
PARAMETER {
	v (mV)
	dt (ms)
	gkbar=.036 (S/cm2)
	ek = -77 (mV)
	celsius = 16 (degC)
}
STATE {
	n
}
ASSIGNED {
	ik (mA/cm2)
	inf
}
LOCAL	fac

INITIAL {
	rate(v*1(/mV))
	n = inf
}

BREAKPOINT {
	SOLVE states
	ik = gkbar*n*n*n*n*n*n*(v - ek)
}

PROCEDURE states() {	: exact when v held constant
	rate(v*1(/mV))
	n = n + fac*(inf - n)
	VERBATIM
	return 0;
	ENDVERBATIM
}

UNITSOFF
FUNCTION alp(v(mV)) { LOCAL q10
	v = -v - 65
	q10 = 3^((celsius - 6.3)/10)
	alp = q10 * .01*expM1(v + 10, 10)
}

FUNCTION bet(v(mV)) { LOCAL q10
	v = -v - 65
	q10 = 3^((celsius - 6.3)/10)
	bet = q10 * .125*exp(v/80)
}

FUNCTION expM1(x,y) {
        if (fabs(x/y) < 1e-6) {
                expM1 = y*(1 - x/y/2)
        }else{
                expM1 = x/(exp(x/y) - 1)
        }
}


PROCEDURE rate(v) {LOCAL a, b, tau :rest = -70
	TABLE inf, fac DEPEND dt, celsius FROM -100 TO 100 WITH 200
		a = alp(v)  b=bet(v)
		tau = 1/(a + b)
		inf = a/(a + b)
		fac = (1 - exp(-dt/tau))
}
UNITSON



