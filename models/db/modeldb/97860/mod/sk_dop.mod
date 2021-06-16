NEURON {
	SUFFIX sk_dop
	USEION k READ ek WRITE ik
        USEION ca READ cai
        RANGE  gbar,gkahp,ik
        GLOBAL inf,tau
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
	(S) = (siemens)
}

PARAMETER {
	celsius = 6.3	(degC)
	gbar = 1	(S/cm2)
        n = 4
        cai = 50.e-6	(mM)
        a0 = 1.3e13	(1/ms-mM-mM-mM-mM)	:b0/(1.4e-4^4)
        b0 = .5e-2	(1/ms)			:0.5/(0.100e3)
}

STATE {	w }

ASSIGNED {
	ik	(mA/cm2)
        g	(S/cm2)
        inf
        tau	(ms)
	a	(1/ms)
        v	(mV)
        ek	(mV)
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	g = gbar*w
	ik = g*(v-ek)
}

INITIAL {
	rate(cai)
	w=inf
}

DERIVATIVE state {
	rate(cai)
	w' = (inf - w)/tau
}

PROCEDURE rate(cai (mM)) {
	a = a0*cai^4
	tau = 1/(a + b0)
        inf = a*tau
}
