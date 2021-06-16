TITLE Delayed rectifier potassium channel

NEURON {
	SUFFIX KDr
	USEION k READ ek WRITE ik
	RANGE gkdrbar, ik, g, i
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S)  = (siemens)
}

PARAMETER {
	gkdrbar 	= .8 	(mho/cm2)
}

ASSIGNED {
	v 		(mV)
	ek          (mV)
	ik		(mA/cm2)
	i (mA/cm2)
	g (S/cm2)
	ninf ntau
}

STATE {
	n
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gkdrbar * n * n * n * n
	i = g * (v - ek)
	ik = i
}

INITIAL {
	:rates(v)
	:n = ninf
	n = 0.1239
}

DERIVATIVE states {
	rates(v)
	n' = (ninf - n) / ntau
}

PROCEDURE rates(v(mV)) {
	ntau = 5 / (exp((v + 50)/40) + exp(-(v + 50)/50))
	ninf = 1 / ((exp(-(v + 31)/15)) + 1)
}




