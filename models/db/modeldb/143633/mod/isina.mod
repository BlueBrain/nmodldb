TITLE initial segment Na channel

: M. Birdno added vtrap functions on 1/30/08
:

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX isina
	NONSPECIFIC_CURRENT ina
	RANGE gnabar, ena
	RANGE m_inf, h_inf
	RANGE tau_m, tau_h
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gnabar	= 0.4 	(mho/cm2)
	ena	= 45	(mV)
	celsius		(degC)
	dt              (ms)
	v               (mV)
	vtraub	= -70	(mV)
}

STATE {
	m h
}

ASSIGNED {
	ina	(mA/cm2)
	m_inf
	h_inf
	tau_m (ms)
	tau_h (ms)
	tadj
}


BREAKPOINT {
	SOLVE states METHOD euler
	ina = gnabar * m*m*m*h * (v - ena)
}


DERIVATIVE states {   : exact Hodgkin-Huxley equations
	evaluate_fct(v)
	m' = (m_inf - m) / tau_m
	h' = (h_inf - h) / tau_h
}

UNITSOFF
INITIAL {
	tadj = 3.0 ^ ((celsius-36)/ 10 )
	evaluate_fct(v)
	m = m_inf
	h = h_inf
}

PROCEDURE evaluate_fct(v(mV)) { LOCAL a,b,v2

	v2 = v - vtraub : convert to traub convention

	a = vtrap1(v2)
	b = vtrap2(v2)
	tau_m = 1 / (a + b) / tadj
	m_inf = a / (a + b)

	a = 0.128 * Exp((17-v2)/18)
	b = 4 / ( 1 + Exp((40-v2)/5) )
	tau_h = 1 / (a + b) / tadj
	h_inf = a / (a + b)
}

UNITSON

FUNCTION vtrap1(x) {
	if (fabs(13-x) < 1e-6) {
		vtrap1 = 0.32*4
	}else{
		vtrap1 =  0.32 * (13-x) / ( Exp((13-x)/4) - 1)
	}
}

FUNCTION vtrap2(x) {
	if (fabs(x-40) < 1e-6) {
		vtrap2 = 0.28*5
	}else{
		vtrap2 =  0.28 * (x-40) / ( Exp((x-40)/5) - 1)
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

