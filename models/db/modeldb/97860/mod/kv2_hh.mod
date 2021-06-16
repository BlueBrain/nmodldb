: HH-style Kv2.1 channel model w/ inactivation
:
: Fits by DP Mohapatra and Josh Held
: 1/7/2005 update! - Adjusted vhn and vcn - Josh Held

NEURON {
	SUFFIX kv2_hh
	USEION k READ ek WRITE ik
	RANGE g, ninf, tn, hinf, th, ik, gbar
	GLOBAL vhn, vcn, vhh, vch, p
	GLOBAL Ctn, vhtn, atn, btn, tn0
	GLOBAL Cth, vhth, ath, bth, th0
}

UNITS {
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER {
	gbar	= 1	(S/cm2)

	:vhn	= 18	(mV)
	vhn	= 17.5	(mV)
	:vcn	= -18	(mV)
	vcn	= -10	(mV)

	vhh	= -25	(mV)
	vch	= -12	(mV)
	p	= .26	

	:Ctn	= 80	(ms)
	:vhtn	= -10	(mV)
	:atn	= 14	(mV)
	:btn	= 20	(mV)
	:tn0	= 5	(ms)

	Ctn	= 5	(ms)
	vhtn	= -30	(mV)
	atn	= 14	(mV)
	btn	= 20	(mV)
	tn0	= 5	(ms)

	
	Cth	= 500	(ms)
	vhth	= 50	(mV)
	ath	= 20	(mV)
	bth	= 20	(mV)
	th0	= 800	(ms)

	Cq10 	= 4
	celsius		(degC)
}

ASSIGNED {
	g       (S/cm2)
	v	(mV)
	ninf
	hinf
	tn	(ms)
	th	(ms)
	ik	(mA/cm2)
	ek	(mV)
}

STATE {
	n
	h
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gbar*n*h
	ik = g*(v-ek)
}

DERIVATIVE states{
	values()
	n' = (ninf - n)/tn
	h' = (hinf - h)/th
}

INITIAL {
	values()
	n = ninf
	h = hinf
}

PROCEDURE values() {LOCAL q10
	q10 = Cq10^((celsius-23 (degC))/10 (degC))
	ninf = 1/(1 + exp((v - vhn)/vcn))
	hinf = (1-p)/(1 + exp(-(v - vhh)/vch)) + p
	tn = q10*Ctn/(exp((v-vhtn)/atn) + exp(-(v-vhtn)/btn)) + tn0
	th = q10*Cth/(exp((v-vhth)/ath) + exp(-(v-vhth)/bth)) + th0
}