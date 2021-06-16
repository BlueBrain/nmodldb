TITLE Id and Ir currents of the Huber-Braun Model
: Braun et al. Int J Bifurcation and Chaos 8(5):881-889 (1998)
: Fast Na+ and K+ currents responsible for action potentials
: 
: Written by Patricio Orio, Jul 2006
:

NEURON {
	SUFFIX dr
	USEION na READ ena WRITE ina
	USEION k READ ek WRITE ik
	RANGE gd, gr
	RANGE td, tr
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gd  = 0.0025  (mho/cm2)
	gr  = 0.0028  (mho/cm2)
	V0d = -25      (mV)
    V0r = -25       (mV)
    sd = 0.25       (/mV)
    sr = 0.25       (/mV)
    tr = 1.5         (ms)
}

STATE {
	ad ar
}

ASSIGNED {
	celsius	(degC)
	ina     (mA/cm2)
	ik      (mA/cm2)
    v       (mV)
	ena     (mV)
    ek      (mV)
    rho     (1)
    arinf
}

INITIAL {
    rho = 1.3^((celsius - 25 (degC))/10(degC))
    ar = 1/(1+exp(-sr*(v - V0r)))
    ad = 1/(1+exp(-sd*(v - V0d)))
}

BREAKPOINT {
	SOLVE states METHOD cnexp
    rho = 1.3^((celsius - 25 (degC))/10(degC))
	ad = 1/(1+exp(-sd*(v - V0d)))
    ina = rho * gd * ad * (v - ena)
	ik  = rho * gr * ar * (v - ek)
}

DERIVATIVE states {
    LOCAL phi
    phi = 3^((celsius - 25 (degC))/ 10 (degC))
    arinf = 1/(1+exp(-sr*(v - V0r)))
    ar' = phi*(arinf - ar)/tr
}

