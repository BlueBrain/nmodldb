TITLE slow K current

COMMENT
12/1/2005 NTC Made compatible with adaptive integration
Unused stuff removed
ENDCOMMENT

: modified by Steven Prescott based on current described below
: Prescott and De Koninck. 2005. J Neurosci 25: 4743-4754
: Slow potassium current in spinal lamina I neurons
: produces slow component of biphasic AHP between spikes
:
: original current described below...
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

NEURON {
	SUFFIX Ks
	USEION k READ ek WRITE ik
	RANGE gkbar, vtraub, rate_change
	RANGE n_inf
	RANGE tau_n
	RANGE ik 
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gkbar		= .0003 	(mho/cm2)
	ek				(mV)
	celsius			(degC)
	v               		(mV)
	vtraub	= -55		(mV)	: adjusts threshold
	rate_change = 0.1			: slows kinetics when <1
}

STATE {
	n
}

ASSIGNED {
	ik	(mA/cm2)
	n_inf
	tau_n (ms)
	tadj
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	ik  = gkbar * n * (v - ek)
}

DERIVATIVE states {
	evaluate_fct(v)
	n' = (n_inf-n)/tau_n
}

UNITSOFF
INITIAL {
:
:  Q10 was assumed to be 3
:
	tadj = 3.0 ^ ((celsius-36)/ 10 )
	evaluate_fct(v)
	n= n_inf
}

PROCEDURE evaluate_fct(v(mV)) { LOCAL a,b,v2

	v2 = v - vtraub : convert to traub convention

	a = 0.032 * (15-v2) / ( exp((15-v2)/5) - 1)
	b = 0.5 * exp((10-v2)/40)
	tau_n = 1 / (a + b) / (tadj*rate_change)
	n_inf = a / (a + b)
}

UNITSON
