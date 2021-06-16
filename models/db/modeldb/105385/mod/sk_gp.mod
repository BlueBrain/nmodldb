NEURON {
	SUFFIX sk_gp
	USEION k READ ek WRITE ik
        USEION ca READ cai
        RANGE  gbar,gkahp,ik
        GLOBAL inf,tau
	GLOBAL Cq10
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
	(S) = (siemens)
}

PARAMETER {
	gbar = 1	(S/cm2)
        n = 4
        cai = 50.e-6	(mM)
        a0 = 1.3e13	(1/ms-mM-mM-mM-mM)	:b0/(1.4e-4^4)
        b0 = .5e-2	(1/ms)			:0.5/(0.100e3)
	        celsius (degC)
	Cq10 = 3
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
	LOCAL q10:, cai_eff
	:if (cai<1e-12) {
	:	cai_eff = 1e-12
	:} else {
	:	cai_eff = cai
	:}
	q10 = Cq10^((celsius - 22 (degC))/10 (degC) )
	a = a0*cai^4
	:a = a0*cai_eff^4
	tau = q10/(a + b0)
	inf = a/(a + b0)
}
