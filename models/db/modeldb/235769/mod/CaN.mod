TITLE N-type Calcium channel

NEURON {
	SUFFIX CaN
	USEION ca READ cai, cao WRITE ica
	RANGE gcanbar, ica, g, i, eca
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S)  = (siemens)
}

PARAMETER {
	gcanbar  = .01 (mho/cm2)
	mtau	= 15	(ms)
	htau	= 50	(ms)
	R = 8.31441 (VC/Mol/K)
	T = 309.15 (k) 
	Z = 2
	F = 96485.309 (/C)
}

ASSIGNED {
	v 		(mV)
	eca 	(mv)
	ica		(mA/cm2)
	i 		(mA/cm2)
	g 		(S/cm2)
	minf hinf
	cai cao	(mM)
}

STATE {
	m h
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	eca = ((1000*R*T)/(Z*F))*log(cao/cai)
	g = gcanbar * m * m * h
	i = g * (v - (eca-70))
	ica = i
}

INITIAL {
	m = 0.004199
	h = 0.9219
}

DERIVATIVE states {
	rates(v)
	m' = (minf - m) / mtau
	h' = (hinf - h) / htau
}

PROCEDURE rates(v(mV)) {
	minf = 1 / (1 + exp((v+25)/-5))
	hinf = 1 / (1 + exp((v+43)/5))
}



