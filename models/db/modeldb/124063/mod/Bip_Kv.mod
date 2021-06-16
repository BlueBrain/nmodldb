
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX IKv
	USEION k READ ek WRITE ik
	RANGE gkvbar
	RANGE m_inf, h_inf
	RANGE tau_m, tau_h
	RANGE m_exp, h_exp

}


UNITS {
	(molar) = (1/liter)
	(mM) = (millimolar)
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {
	gkvbar	= 0.000318	(mho/cm2)
	ek		(mV)
	dt              (ms)
	v               (mV)

}

STATE {
	m h 
}

INITIAL {
: The initial values were determined at a resting value of -66.3232 mV in a single-compartment
: at -60 mV
        m = 0.0345
        h = 0.8594
}

ASSIGNED {
	 ik	(mA/cm2)
       m_inf 
	 h_inf
	 tau_m 
	 tau_h
	 m_exp 
       h_exp

}

BREAKPOINT {
	SOLVE states
	ik = gkvbar * m*m*m *h* (v - ek)

}

PROCEDURE states() {	: exact when v held constant
	evaluate_fct(v)
	m = m + m_exp * (m_inf - m)
	h = h + h_exp * (h_inf - h)

	VERBATIM
	return 0;
	ENDVERBATIM

}

UNITSOFF

PROCEDURE evaluate_fct(v(mV)) { LOCAL am,bm,ah,bh
	
:IA channel

 am = (400) / ((exp(-(v-15)/36)) + 1)
 bm = (exp((-v/13)))
 tau_m = 1 / (am + bm)
 m_inf = am * tau_m
 m_exp = 1 - exp(-dt/tau_m)

 ah = 0.0003 * (exp((-v/7)))
 bh = ((80) / ((exp(-(v+115)/15)) + 1))+0.02
 tau_h = 1 / (ah + bh)
 h_inf = ah * tau_h
 h_exp = 1 - exp(-dt/tau_h)

}

UNITSON
