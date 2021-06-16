TITLE ICa for horizontal cell
: Calcium current (ICa) for horizontal cells
: 
: Based on parameters of Aoyama et al. (2000)


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX HzICa
	USEION ca WRITE ica
	RANGE gbar
	RANGE m_inf
	RANGE tau_m
	RANGE m_exp
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
    gbar = 8.4906e-05 (mho/cm2) : 9 ns total
    eca  = 53.040  (mV)         : 2 mM outside, 30 uM inside
}

STATE {
    m
}

ASSIGNED {
    v       (mV)
	ica     (mA/cm2)
    celsius (degC)
    dt      (ms)

    m_inf
	tau_m
	m_exp
	tadj
}

BREAKPOINT {
	SOLVE states
	ica = gbar * m*m*m*m * (v - eca)
}

PROCEDURE states() {
    : exact when v held constant
	evaluate_fct(v)
	m = m + m_exp * (m_inf - m)
}

UNITSOFF

INITIAL {
	m = 0.059
    tadj = 3.0 ^ ((celsius-25)/10)  : correction for physio temp
}

PROCEDURE evaluate_fct(v(mV)) { LOCAL a, b, v2
    a = (240 * (68-v)) / ( exp( (68-v)/21 ) - 1)
    b = 800 / ( exp( (55+v)/55 ) + 1 )
    tau_m = 1 / (a + b) / tadj
	m_inf = a / (a + b)

	m_exp = 1 - exp(-dt/tau_m)
}

UNITSON

