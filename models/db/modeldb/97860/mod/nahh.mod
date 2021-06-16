: Na HH model

NEURON {
	SUFFIX nahh
	USEION na READ ena WRITE ina
	RANGE i, g, gbar, minf, hinf, tm, th
    GLOBAL vhm, vcm, vhh, vch, p
    GLOBAL Ctm, vhtm, atm, btm, tm0
    GLOBAL Cth, vhth, ath, bth, th0
}

UNITS {
	(S)     = (siemens)
	(mV)    = (millivolt)
	(mA)    = (milliamp)
}

PARAMETER {
	gbar = 1	(S/cm2)

	vhm = -25	(mV)
	vcm = 7		(mV)
	vhh = -52	(mV)
	vch = 5		(mV)
	p = 0
	Ctm = 0		(ms)
	vhtm = -40	(mV)
	atm	= 10	(mV)
	btm = 10	(mV)
	tm0 = 0.3	(ms)
	Cth = 0	(ms)
	vhth = -60	(mV)
	ath	= 10	(mV)
	bth = 10	(mV)
	th0 = 5	(ms)
}

ASSIGNED {
	v   	(mV)
	minf
	hinf
	tm  	(ms)
	th		(ms)
	ina		(mA/cm2)
    g		(S/cm2)
    ena		(mV)
    celsius	(degC)
}

STATE {
	m
	h
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gbar*(m^3)*h
    ina = g*(v-ena)
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
	q10 = 3^((celsius-23 (degC))/10 (degC))
	minf = 1/(1+exp(-(v-vhm)/vcm))
	hinf = p + ((1-p)/(1+exp((v-vhh)/vch)))
	tm = q10*((Ctm/(exp((v-vhtm)/atm) + exp(-(v-vhtm)/btm))) + tm0)
	th = q10*((Cth/(exp((v-vhth)/ath) + exp(-(v-vhth)/bth))) + th0)
}