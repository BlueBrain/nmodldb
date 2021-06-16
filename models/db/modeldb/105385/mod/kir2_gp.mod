NEURON {
	SUFFIX kir2_gp
	USEION k READ ek WRITE ik
	RANGE g, ninf, ik, gbar
	GLOBAL vh, vc
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
}

PARAMETER {
	gbar = 1	(S/cm2)
	ek		(mV)
	vh = -90	(mV)
	vc = 12.1	(mV)
}

ASSIGNED {
	v		(mV)
	ninf
	ik		(mA/cm2)		
	g		(S/cm2)
}

STATE {	
}

BREAKPOINT {
	values()
	g = gbar*ninf
	ik = g*(v-ek)	
}

INITIAL {
	values()
}

PROCEDURE values() {
	ninf = 1/(1 + exp((v - vh)/vc))
}