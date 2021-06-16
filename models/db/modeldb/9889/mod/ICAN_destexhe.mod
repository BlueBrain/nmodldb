: $Id: ICAN_destexhe.mod,v 1.11 1994/10/27 21:27:31 billl Exp $
TITLE Slow Ca-dependent cation current
:
:   Ca++ dependent nonspecific cation current ICAN
:   Differential equations
:
:   Model of Destexhe, 1992.  Based on a first order kinetic scheme
:      <closed> + n cai <-> <open>	(alpha,beta)
:
:   Following this model, the activation fct will be half-activated at 
:   a concentration of Cai = (beta/alpha)^(1/n) = cac (parameter)
:   The mod file is here written for the case n=2 (2 binding sites)
:   ---------------------------------------------
:
:   Kinetics based on: Partridge & Swandulla, TINS 11: 69-72, 1988.
:
:   This current has the following properties:
:      - inward current (non specific for cations Na, K, Ca, ...)
:      - activated by intracellular calcium
:      - NOT voltage dependent
:
:   Written by Alain Destexhe, Salk Institute, Dec 7, 1992
:

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX ican
	USEION other WRITE iother VALENCE 1
	USEION ca READ cai
        RANGE gbar, i
	GLOBAL m_inf, tau_m, beta, cac, taumin, erev
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
}


PARAMETER {
	v		(mV)
	celsius	= 36	(degC)
	erev = -20	(mV)
	cai 	= .00005	(mM)	: initial [Ca]i = 50 nM
	gbar	= 1e-5	(mho/cm2)
	beta	= 2.5	(1/ms)		: backward rate constant
	cac	= 1e-4	(mM)		: middle point of activation fct
	taumin	= 0.1	(ms)		: minimal value of time constant
}


STATE {
	m
}

INITIAL {
	evaluate_fct(v,cai)
	m = m_inf
:
:  activation kinetics are assumed to be at 22 deg. C
:  Q10 is assumed to be 3
:
	tadj = 3.0 ^ ((celsius-22.0)/10)

}

ASSIGNED {
	i	(mA/cm2)
	iother	(mA/cm2)
	m_inf
	tau_m	(ms)
	tadj
}

BREAKPOINT { 
	SOLVE states METHOD runge
	i = gbar * m*m * (v - erev)
	iother = i
}

DERIVATIVE states { 
	evaluate_fct(v,cai)

	m' = (m_inf - m) / tau_m
}

UNITSOFF

PROCEDURE evaluate_fct(v(mV),cai(mM)) {  LOCAL alpha2

	alpha2 = beta * (cai/cac)^2

	tau_m = 1 / (alpha2 + beta) / tadj
	m_inf = alpha2 / (alpha2 + beta)

        if(tau_m < taumin) { tau_m = taumin } 	: min value of time cst
}
UNITSON
