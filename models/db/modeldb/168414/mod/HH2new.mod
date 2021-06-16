TITLE Hippocampal HH channels

COMMENT
12/1/2005 NTC Made compatible with adaptive integration
Unused stuff removed
ENDCOMMENT

: Fast Na+ and K+ currents responsible for action potentials
: Iterative equations
:
: Equations modified by Traub, for Hippocampal Pyramidal cells, in:
: Traub & Miles, Neuronal Networks of the Hippocampus, Cambridge, 1991
:
: range variable vtraub adjust threshold
:
: Written by Alain Destexhe, Salk Institute, Aug 1992
:
: Modifications by Arthur Houweling for use in MyFirstNEURON
:
: Addition of Vsm by Steven Prescott

NEURON {
	SUFFIX HH2new
	USEION na READ ena WRITE ina
	USEION k READ ek WRITE ik
	RANGE gnabar, gkbar, vtraub, vsm
	RANGE m_inf, h_inf, n_inf
	RANGE tau_m, tau_h, tau_n
	RANGE ik, ina 
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gnabar	= .1 	(mho/cm2)
	gkbar		= .06 (mho/cm2)

	ena		(mV)
	ek		(mV)
	celsius		(degC)
	v               (mV)
	vtraub	= -55	(mV)	: adjusts threshold
	vsm		= -5 (mV)	: collapses activation curve as increasingly -ve
}

STATE {
	m h n
}

ASSIGNED {
	ina	(mA/cm2)
	ik	(mA/cm2)
	il	(mA/cm2)
	m_inf
	h_inf
	n_inf
	tau_m (ms)
	tau_h (ms)
	tau_n (ms)
	tadj
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	ina = gnabar * m*m*m*h * (v - ena)
	ik  = gkbar * n*n*n*n * (v - ek)
}


DERIVATIVE states {   : use this for exact Hodgkin-Huxley equations
	evaluate_fct(v)
	m' = (m_inf - m) / tau_m
	h' = (h_inf - h) / tau_h
	n' = (n_inf - n) / tau_n
}

UNITSOFF
INITIAL {
:
:  Q10 was assumed to be 3 for both currents
:
	tadj = 3.0 ^ ((celsius-36)/ 10 )
	evaluate_fct(v)
	m= m_inf
	h= h_inf
	n= n_inf
}

PROCEDURE evaluate_fct(v(mV)) { LOCAL a,b,v2

	v2 = v - vtraub : convert to traub convention

	a = 0.32 * (vsm+13-v2) / ( exp((vsm+13-v2)/4) - 1)
	b = 0.28 * (vsm+v2-40) / ( exp((vsm+v2-40)/5) - 1)
	tau_m = 1 / (a + b) / tadj
	m_inf = a / (a + b)

	a = 0.128 * exp((17-v2)/18)
	b = 4 / ( 1 + exp((40-v2)/5) )
	tau_h = 1 / (a + b) / tadj
	h_inf = a / (a + b)

	a = 0.032 * (15-v2) / ( exp((15-v2)/5) - 1)
	b = 0.5 * exp((10-v2)/40)
	tau_n = 1 / (a + b) / tadj
	n_inf = a / (a + b)
}

UNITSON
