: $Id: ICAN_voltdep.mod,v 1.4 1994/04/14 01:28:04 billl Exp $
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
:      - voltage-dependent: a voltage-dependence of ICAN was described 
:        for some cells (cfr. Partridge & Swandulla).  In nRt cells,
:        the study of Bal & McCormick strongly suggests that ICAN 
:        decreases with hyperpolarization.
:
:   The voltage-dependence of ICAN is assumed to be monoexponential
:   with voltage for the two rate constants alpha and beta, such as
:   m_inf is a sigmoid fct which becomes null with hyperpolarization.
:   So ICAN, is a noninactivating current, activated by Ca++ and 
:   depolarization...
:        
:
:   Written by Alain Destexhe, Salk Institute, Dec 7, 1992
:

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX icanv
	USEION n READ en WRITE in VALENCE 1
	USEION ca READ cai
        RANGE gbar
	GLOBAL 	m_inf, tau_m, cac, taumin, vact, vtau
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
	en		(mV)
	cai 	= .00005	(mM)	: initial [Ca]i = 50 nM
	gbar	= 1e-5	(mho/cm2)
	cac	= 1e-4	(mM)		: middle point of activation fct
	taumin	= 0.1	(ms)		: minimal value of time constant
	vact	= -64	(mV)		: half-activation voltage for activ
	vtau	= -92	(mV)		: voltage for time cst exponential
}


STATE {
	m
}

INITIAL {
	evaluate_fct(v,cai)
	m = m_inf
}


ASSIGNED {
	in	(mA/cm2)
	m_inf
	tau_m	(ms)
}

BREAKPOINT { 
	SOLVE states
	in = gbar * m*m * (v - en)
}

DERIVATIVE states { 
	evaluate_fct(v,cai)

	m' = (m_inf - m) / tau_m
}

UNITSOFF
PROCEDURE evaluate_fct(v(mV),cai(mM)) {  LOCAL cc,tadj
:
:  activation kinetics are assumed to be at 22 deg. C
:  Q10 is assumed to be 3
:
:
:
	tadj = 3 ^ ((celsius-22.0)/10)

	cc = (cai/cac)^2

	m_inf = 1 / (1 + exp(-(v-vact)/2) / cc )

	tau_m = exp((v-vtau)/4) / (1 + cc*exp((v-vact)/2) ) / tadj

        if(tau_m < taumin) { tau_m = taumin } 	: min value of time cst
}
UNITSON
