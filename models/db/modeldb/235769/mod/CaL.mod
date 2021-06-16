TITLE L-type Calcium channel

NEURON {
	POINT_PROCESS CaL
	USEION ca WRITE ica
	RANGE gcalbar, g, eca, i, ica
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	gcalbar  = 7e-07 (uS)
	mtau	 = 60	    (ms)
	eca      = 60
}

ASSIGNED {
	ica		(nA)
	v		(mV)
	g		(uS)
	i		(nA)
	minf
}

STATE {
	m
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gcalbar * m
	i =  g * (v - eca)
	ica=i
}

INITIAL {
	rates(v)
	m = minf
}

DERIVATIVE states {
	rates(v)
	m' = (minf - m) / mtau
}

PROCEDURE rates(v(mV)) {
	minf = 1 / (1 + exp(-(v + 43)/6))
}