TITLE Sodium channels

COMMENT
-----------------------------------------------------------------------------
	Na current for action potentials
	--------------------------------

  - fast sodium current
  - iterative equations

  Model of INa for hippocampal pyramidal cells, from
  Traub & Miles, Neuronal Networks of the Hippocampus, Cambridge, 1991

  Added a shift parameter for inactivation


  Written by Alain Destexhe, Laval University, 1996

  modified by Philipp Vetter 28.7.1998

  inaT <-> na
  gnabar <-> gbar

  TABLEs in order to avoid division by zero, Arnd Roth 12.9.99
-----------------------------------------------------------------------------
ENDCOMMENT


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX na
	USEION na READ ena WRITE ina
	RANGE gbar, vtraub, shift
	RANGE m_inf, h_inf
	RANGE tau_m, tau_h
	RANGE m_exp, h_exp
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gbar	 	(pS/um2)	: max conductance
	vtraub	= -63	(mV)		: adjusts threshold
	shift	= 0	(mV)		: inactivation shift
	ena	= 50	(mV)
	celsius = 36    (degC)
	dt              (ms)
	v               (mV)
	vmin = -120	(mV)
	vmax = 100	(mV)
}

STATE {
	m h
}

ASSIGNED {
	ina	(mA/cm2)
	m_inf
	h_inf
	tau_m	(ms)
	tau_h	(ms)
	m_exp
	h_exp
	tadj
}


BREAKPOINT {
	SOLVE states
	ina = (1e-4)*gbar * m*m*m*h * (v - ena)
}


:DERIVATIVE states {
:	evaluate_fct(v)
:	m' = (m_inf - m) / tau_m
:	h' = (h_inf - h) / tau_h
:}

PROCEDURE states() {	: exact when v held constant
	evaluate_fct(v)
	m = m + m_exp * (m_inf - m)
	h = h + h_exp * (h_inf - h)
}

UNITSOFF
INITIAL {
	m = 0
	h = 0
:
:  Q10 was assumed to be 2.3 for both currents
:
:  original measurements at room temperature

	tadj = 2.3 ^ ((celsius-23)/ 10 )
	evaluate_fct(v)
	m = m_inf
	h = h_inf
}



FUNCTION myexp(x) {
	if (x < -100) {
	myexp = 0
	}else{
	myexp = exp(x)
	}
}

PROCEDURE evaluate_fct(v(mV)) { LOCAL a,b,v2

        TABLE m_inf, m_exp, h_inf, h_exp
	DEPEND dt, celsius, shift, vtraub
	
	FROM vmin TO vmax WITH 199

	v2 = v - vtraub 	: convert to traub convention

	a = 0.32 * (13-v2) / ( myexp((13-v2)/4) - 1)
	b = 0.28 * (v2-40) / ( myexp((v2-40)/5) - 1)
	tau_m = 1 / (a + b)
	m_inf = a / (a + b)

	v2 = v2 - shift		: inactivation shift
				: about -10 mV for neocortical pyr cells

	a = 0.128 * myexp((17-v2)/18)
	b = 4 / ( 1 + myexp((40-v2)/5) )
	tau_h = 1 / (a + b)
	h_inf = a / (a + b)

	m_exp = 1 - myexp(-dt/tau_m)
	h_exp = 1 - myexp(-dt/tau_h)
}

UNITSON
