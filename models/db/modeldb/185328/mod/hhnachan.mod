TITLE HH sodium channel
: Hodgkin - Huxley squid sodium channel

NEURON {
	SUFFIX HHna
	USEION na READ ena WRITE ina
	RANGE gnabar, ina
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
	celsius = 6.3	(degC)
	dt (ms)
	gnabar=.120 (S/cm2)
	ena = 50 (mV)
}
STATE {
	m h
}
ASSIGNED {
	ina (mA/cm2)
	inf[2]
}
LOCAL	fac[2]

INITIAL {
	rate(v*1(/mV))
	m = inf[0]
	h = inf[1]
}

BREAKPOINT {
	SOLVE states
	ina = gnabar*m*m*m*h*(v - ena)
}

PROCEDURE states() {	: exact when v held constant
	rate(v*1(/mV))
	m = m + fac[0]*(inf[0] - m)
	h = h + fac[1]*(inf[1] - h)
	VERBATIM
	return 0;
	ENDVERBATIM
}

UNITSOFF
FUNCTION alp(v(mV),i) { LOCAL a,b,c,q10 :rest = -70  order m,h
	v = -v - 65  :convert to hh convention
	q10 = 3^((celsius - 6.3)/10)
	if (i==0) {
		alp = q10*.1*expM1(v + 25, 10)
	}else if (i==1){
		alp = q10*.07*exp(v/20)
	}
}

FUNCTION bet(v,i) { LOCAL a,b,c,q10 :rest = -70  order m,h
	v = -v - 65
	q10 = 3^((celsius - 6.3)/10)
	if (i==0) {
		bet = q10* 4*exp(v/18)
	}else if (i==1){
		bet = q10*1/(exp(.1*v + 3) + 1)
	}
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
	FROM i=0 TO 1 {
		a = alp(v,i)  b=bet(v,i)
		tau = 1/(a + b)
		inf[i] = a/(a + b)
		fac[i] = (1 - exp(-dt/tau))
	}
}
UNITSON


