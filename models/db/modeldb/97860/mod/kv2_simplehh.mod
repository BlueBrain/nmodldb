: HH-style non-inactivating Kv2 channel model

NEURON {
	SUFFIX kv2_simplehh
	USEION k READ ek WRITE ik
	RANGE g, ninf, tn, ik, gbar
	GLOBAL vhn, vcn
	GLOBAL Ctn, vhtn, atn, btn, tn0
}

UNITS {
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER {
	gbar	= 1	(S/cm2)
	
	vhn	= 0	(mV)
	vcn	= -8	(mV)	

	Ctn	= 5	(ms)
	vhtn	= -30	(mV)
	atn	= 14	(mV)
	btn	= 20	(mV)
	tn0	= 5	(ms)

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
	tn = q10*Ctn/(exp((v-vhtn)/atn) + exp(-(v-vhtn)/btn)) + tn0
}