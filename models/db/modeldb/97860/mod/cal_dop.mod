NEURON {
	SUFFIX cal_dop
	USEION ca READ cai,cao WRITE ica
        RANGE  gbar,ica,g, vf
	GLOBAL vhm, vcm, voff
	GLOBAL Ctm, atm, btm, tm0, vhtm
        GLOBAL minf,tau
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	FARADAY = (faraday)  (kilocoulombs)
	R = (k-mole) (joule/degC)
	KTOMV = .0853 (mV/degC)
	(S) = (siemens)
	(mM) = (milli/liter)
}

PARAMETER {
	v		(mV)
	voff = 0	(mV)
	celsius	(degC)
	gbar = .003 	(S/cm2)
	ki = .001 	(mM)
	cai 		(mM)
	cao 		(mM)
        tfa = 1
	vhm = -35	(mV)
	vcm = 8	(mV)
	Ctm = 3		(ms)
	atm = 12	(mV)
	btm = 11	(mV)
	tm0 = 0		(ms)
	vhtm = -2	(mV)
}

STATE { m }

ASSIGNED {
	vf	(mV)
	ica	(mA/cm2)
        g	(S/cm2)
        minf
        tau	(ms)
	a	(1/ms)
}

BREAKPOINT {
	vf = v-voff
	SOLVE state METHOD cnexp
	g = gbar*m*m*h2(cai)
	ica  = g*ghk(vf,cai,cao)
}

INITIAL {
	vf = v-voff
	rate(vf)
	m = minf
}

FUNCTION h2(cai(mM)) {
	h2 = ki/(ki+cai)
}

FUNCTION ghk(vf(mV), ci(mM), co(mM)) (mV) {
        LOCAL nu,f

        f = KTF(celsius)/2
        nu = vf/f
        ghk=-f*(1. - (ci/co)*exp(nu))*efun(nu)
}

FUNCTION KTF(celsius (degC)) (mV) {
        KTF = ((25 (mV) /293.15 (degC) )*(celsius + 273.15 (degC) ))
}

FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}

FUNCTION alp(vf(mV)) (1/ms) {
	TABLE FROM -150 TO 150 WITH 200
	alp = 15.69 (1/mV-ms) *(-1.0*vf+81.5 (mV) )/(exp((-1.0*vf+81.5 (mV) )/10.0 (mV) )-1.0)
}

FUNCTION bet(vf(mV)) (1/ms) {
	TABLE FROM -150 TO 150 WITH 200
	bet = 0.29 (1/ms) *exp(-vf/10.86 (mV) )
}

DERIVATIVE state {
        rate(vf)
        m' = (minf - m)/tau
}

PROCEDURE rate(vf (mV)) { :callable from hoc
        a    = alp(vf)
        :tau  = 1/(tfa*(a + bet(vf)))
        :minf = tfa*a*tau
	tau = Ctm/(exp((vf-vhtm)/atm) + exp((vhtm-vf)/btm)) + tm0
	minf = 1/(1+exp(-(vf-vhm)/vcm))
}
