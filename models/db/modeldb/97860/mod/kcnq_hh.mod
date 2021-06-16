: HH-style KCNQ channel model

NEURON {
	SUFFIX kcnq_hh
	USEION k READ ek WRITE ik
	RANGE g, ninf, tn, ik, gbar
	GLOBAL vhn, vcn
	GLOBAL Ctn, vhtn, vctn, tn0
}

UNITS {
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER {
	gbar	= 1	(S/cm2)
	
	vhn	= -49	(mV)
	vcn	= -9.7	(mV)	

	Ctn	= 100	(ms)
	vhtn = -49	(mV)
	vctn = -9.7	(mV)
	tn0	= 55	(ms)

	Cq10 	= 4
	celsius		(degC)
}

ASSIGNED {
	g       (S/cm2)
	v	(mV)
	ninf
	tn	(ms)
	ik	(mA/cm2)
	ek	(mV)
}

STATE {
	n
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gbar*n
	ik = g*(v-ek)
}

DERIVATIVE states{
	values()
	n' = (ninf - n)/tn
}

INITIAL {
	values()
	n = ninf
}

PROCEDURE values() {LOCAL q10
	q10 = Cq10^((celsius-23 (degC))/10 (degC))
	ninf = 1/(1 + exp((v - vhn)/vcn))
	tn = ((Ctn - tn0)/(1 + exp((v - vhtn)/vctn))) + tn0
}