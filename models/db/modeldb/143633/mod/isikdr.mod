TITLE Delayed Rectifier K channel

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX isikdr
	USEION k READ ek WRITE ik
	RANGE gkdrbar, m_inf, tau_m, i
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gkdrbar	= .01 	(mho/cm2)

:	ek	= -95	(mV)
	celsius		(degC)
	dt              (ms)
	v               (mV)
	vtraub	= -60	(mV)
}

STATE {
	n
}

ASSIGNED {
	ek (mV)
	i	(mA/cm2)
	ik	(mA/cm2)
	n_inf
	tau_n (ms)
	tadj
}


BREAKPOINT {
	SOLVE states METHOD euler
	i = gkdrbar * n*n*n*n * (v - ek)
	ik  = i
}


DERIVATIVE states {   : exact Hodgkin-Huxley equations
	evaluate_fct(v)
	n' = (n_inf - n) / tau_n
}

UNITSOFF
INITIAL {
	tadj = 3^((celsius-36)/10)
	evaluate_fct(v)
	n = n_inf
}

PROCEDURE evaluate_fct(v(mV)) { LOCAL a,b,v2

	v2 = v - vtraub : convert to traub convention

	a = vtrap1(v2)
	b = 0.5 * Exp((10-v2)/40)
	tau_n = 1 / (a + b) / tadj
	n_inf = a / (a + b)
}

UNITSON

FUNCTION vtrap1(x) {
	if (fabs(15-x) < 1e-6) {
		vtrap1 = 0.032*5
	}else{
		vtrap1 =  0.032 * (15-x) / ( Exp((15-x)/5) - 1)
	}
}

FUNCTION Exp(x) {
	if (x < -100) {
		:Exp = 0
	}else{
		if (x > 700) {
			Exp = exp(700)
		}else{
			Exp = exp(x)
		}
	}
}
