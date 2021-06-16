TITLE HH k channel channel
: Hodgkin - Huxley k channel


NEURON {
	SUFFIX B_DR
	USEION k READ ek WRITE ik
	RANGE gkbar, ik
	GLOBAL inf
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
PARAMETER {
	v (mV)
	dt (ms)
	gkbar=.036 (mho/cm2) <0,1e9>
	ek = -84 (mV)
	celsius = 6.3 (degC)
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
	ik = gkbar*n*n*n*n*(v - ek)
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
	v = v
	q10 = 3^((celsius - 6.3)/10)
	alp = q10 * .0075*expM1(-v - 30, 10)
}

FUNCTION bet(v(mV)) { LOCAL q10
	v = v
	q10 = 3^((celsius - 6.3)/10)
	bet = q10 * .1*exp((-v - 46)/31)
}

FUNCTION expM1(x,y) {
        if (fabs(x/y) < 1e-6) {
                expM1 = y*(1 - x/y/2)
        }else{
                expM1 = x/(exp(x/y) - 1)
        }
}


PROCEDURE rate(v) {LOCAL a, b, tau :rest = -70
	TABLE inf, fac DEPEND dt, celsius FROM -150 TO 100 WITH 200
		a = alp(v)  b=bet(v)
		tau = 1/(a+b)
		inf = a/(a + b)
		fac = (1 - exp(-dt/tau))
}
UNITSON
