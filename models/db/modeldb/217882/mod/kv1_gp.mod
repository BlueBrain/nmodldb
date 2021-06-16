:   KV1_GP.MOD
:
:   Kv1.2 channel model using HH-type activation/inactivation
:
:   3/2003
:   Josh Held

NEURON {
	SUFFIX kv1_gp
	USEION k READ ek WRITE ik
	RANGE th, tm, ik, hinf, minf, g, gbar
	GLOBAL p
	GLOBAL vhm, vcm
	GLOBAL vhh, vch
	GLOBAL Cth, vhth, ath, bth, th0
	GLOBAL tm0, Ctm, vhtm, vctm
	GLOBAL th90
	GLOBAL Cq10
}

UNITS {
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER {
	gbar = 1	(S/cm2)
	ek		(mV)
	vhm = -27	(mV)
	vhh = -33.477	(mV)
	vcm = -16	(mV)
	vch = 21.5	(mV)
	Cth = 548.67	(ms)
	vhth = -0.956	(mV)
	ath = 29.013	(mV)
	bth= 100	(mV)
	th0 = 779	(ms)
	tm0 = 3.4	(ms)
	Ctm = 89.2	(ms)
	vhtm = -34.3	(mV)
	vctm = 30.1	(mV)
	p = 0.004
	celsius		(degC)
	Cq10 = 3
}

ASSIGNED {
	v	(mV)
	minf
	hinf
	tm	(ms)
	th	(ms)
	th90	(ms)
	ik	(mA/cm2)
	g	(S/cm2)
}

STATE {
	m
	h
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gbar * (m^2) * h 
	ik = g * (v - ek) 
}

DERIVATIVE states{
	values()
	m' = (minf - m)/tm
	h' = (hinf - h)/th
}

INITIAL {
	values()
	m = minf
	h = hinf
}

PROCEDURE values() {LOCAL q10
	q10 = Cq10^((celsius-22 (degC))/10 (degC))
	minf = 1/(1 + exp((v - vhm)/vcm))
	tm = (1/q10)*(tm0 + Ctm*exp(-((v-vhtm)/vctm)^2))

	hinf = (1-p)/(1 + exp((v - vhh)/vch)) + p
	th = (1/q10)*((Cth/(exp((v-vhth)/ath) + exp(-(v-vhth)/bth))) + th0)
	th90 = (1/q10)*((Cth/(exp((-90-vhth)/ath) + exp(-(-90-vhth)/bth))) + th0)
}





