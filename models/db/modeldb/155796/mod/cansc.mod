TITLE Ca-dependent non-specific cation current
: Original model written by Alain Destexhe, Salk Institute, Dec 7, 1992
: Kinetics based on: Partridge & Swandulla, TINS 11: 69-72, 1988.

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX cansc
	USEION other2 WRITE iother2 VALENCE 1
	USEION ca READ cai
        RANGE gbar, i, g, ratc
	GLOBAL m_inf, tau_m, beta, cac, taumin, erev, x
	THREADSAFE m_inf, tau_m, beta, cac, taumin, erev, x
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
}


PARAMETER {
	v		(mV)
	celsius		(degC)
	erev = -38	(mV)
	cai 		(mM)
	gbar	= 4e-4	(mho/cm2)
	beta	= 1e-3	(1/ms)		: backward rate constant
	cac	= 5e-4	(mM)		: middle point of activation fct
	cas	= 2e-5	(mM)		: middle point of activation fct
	taumin	= -0.1	(ms)		: minimal value of time constant
        ratc    = 1e-1
        x       = 2
}


STATE {
	m
}

INITIAL {
:
:  activation kinetics are assumed to be at 22 deg. C
:  Q10 is assumed to be 3
:

	tadj = 3.0 ^ ((celsius-22.0)/10)

	evaluate_fct(v,cai)
	m = m_inf
}

ASSIGNED {
	i	(mA/cm2)
	iother2	(mA/cm2)
	g       (mho/cm2)
	m_inf
	tau_m	(ms)
	tadj
}

BREAKPOINT { 
	SOLVE states METHOD cnexp
	g = gbar * m*m
	i = g * (v - erev)
	iother2 = i
}

DERIVATIVE states { 
	evaluate_fct(v,cai)

	m' = (m_inf - m) / tau_m
}

UNITSOFF

PROCEDURE evaluate_fct(v(mV),cai(mM)) {  LOCAL alpha2:, tcar
  
	alpha2 = ratc/(1+exp((cac-cai)/cas))
 
	tau_m = 1 / (alpha2 + beta) / tadj
	m_inf = alpha2 / (alpha2 + beta)

        if(tau_m < taumin) { tau_m = taumin } 	: min value of time cst
}
UNITSON
