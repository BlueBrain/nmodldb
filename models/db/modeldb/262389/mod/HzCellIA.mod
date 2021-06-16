TITLE IA for horizontal cell
: Transient outward potassium current (IA) for horizontal cells
: 
: Based on parameters of Aoyama et al. (2000)


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX HzIA
	USEION k READ ek WRITE ik
	RANGE gbar
	RANGE m_inf, h_inf
	RANGE tau_m, tau_h
	RANGE m_exp, h_exp
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
    gbar = 1.4151e-04 (mho/cm2)  : 15 nS total
}

STATE {
    m h
}

ASSIGNED {
    v       (mV)
    ek      (mV)
	ik	    (mA/cm2)
    celsius (degC)
    dt      (ms)

    m_inf
	h_inf
	tau_m
	tau_h
	m_exp
	h_exp
	tadj
}

BREAKPOINT {
	SOLVE states
	ik = gbar * m*m*m*h * (v - ek)
}

PROCEDURE states() {
    : exact when v held constant
	evaluate_fct(v)
	m = m + m_exp * (m_inf - m)
	h = h + h_exp * (h_inf - h)
}

UNITSOFF

INITIAL {
	m = 0.030
	h = 0.998
    tadj = 3.0 ^ ((celsius-25)/10)  : correction for physio temp
}

PROCEDURE evaluate_fct(v(mV)) { LOCAL a, b, v2

    a = 2400 / ( 1 + exp( -(v-50)/28 ) )
    b = 80 * exp( -v/36 )
    tau_m = 1 / (a + b) /tadj
	m_inf = a / (a + b)

    a = exp( -v/60 )
    b = 20 / ( exp( -(v+40)/5 ) + 1)
	tau_h = 1 / (a + b) / tadj
	h_inf = a / (a + b)

	m_exp = 1 - exp(-dt/tau_m)
	h_exp = 1 - exp(-dt/tau_h)
}

UNITSON

