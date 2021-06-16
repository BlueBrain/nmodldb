: A conductance

NEURON {
	SUFFIX KA
	USEION k WRITE ik
	RANGE gbar, minf, mtau, hinf, htau, i, g, m, h
	GLOBAL erev
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
} 

PARAMETER {
	erev = -85		(mV)		: effective Ek
	gbar = 0   		(S/cm2)
 	vmin = -100		(mV)		: for look-up table
	vmax = 100		(mV)
} 


ASSIGNED {
	v 		(mV)
	i 		(mA/cm2)
	ik 		(mA/cm2)
	g			(S/cm2)

	malpha	(/ms)		
	mbeta		(/ms)	
	minf
	mtau 		(ms)

	halpha	(/ms)		
	hbeta		(/ms)
	hinf
	htau 		(ms)
}

STATE { m h }

INITIAL { 
	rates(v)
	m = minf
	h = hinf
}

BREAKPOINT {
  SOLVE states METHOD cnexp
	g = gbar*m*h
	ik = g*(v - erev)
	i = ik
} 

DERIVATIVE states {
	rates(v)
	m' = (minf - m)/mtau
	h' = (hinf - h)/htau
}

PROCEDURE rates(v (mV)) {
TABLE minf, mtau, hinf, htau
FROM vmin TO vmax WITH 199
  malpha = -0.05*(v+25)/(exp(-(v+25)/15)-1)
  mbeta = 0.1*(v+15)/(exp((v+15)/8)-1)
  mtau = 1/(malpha + mbeta)
  minf = malpha/(malpha + mbeta)

  halpha = (1.5e-4)/exp((v+13)/15)
  hbeta = 0.06/(exp(-(v+68)/12)+1)
	htau = 1/(halpha + hbeta)
	hinf = halpha/(halpha + hbeta)
}

