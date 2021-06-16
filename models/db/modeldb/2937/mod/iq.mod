TITLE Q current 

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
        (mM) = (milli/liter)

}

PARAMETER {
        dt (ms)
	v (mV)
	erevq=-55   (mV)
	celsius 	(degC)
	gqbar=.003 (mho/cm2)
        vhalf=-93   (mV)
        a0=0.0009      (/ms)
        zeta=5    (1)
        gq=0.4   (1)
}


NEURON {
	SUFFIX q
	NONSPECIFIC_CURRENT Iq
        RANGE Iq,gqbar
        GLOBAL inf,tau
}

STATE {
        q
}

ASSIGNED {
	Iq (mA/cm2)
        inf
        tau
}

INITIAL {
        rate(v)
        q=inf

}

BREAKPOINT {
	SOLVE state METHOD cnexp
	Iq = gqbar*q*(v-erevq)

}

FUNCTION alp(v(mV)) {
  alp = exp( 1.e-3*zeta*(v-vhalf)*9.648e4/(8.315*(273.16+celsius)))
}

FUNCTION bet(v(mV)) {
  bet = exp(1.e-3*zeta*gq*(v-vhalf)*9.648e4/(8.315*(273.16+celsius))) 
}

DERIVATIVE state { 
        rate(v)
        q' = (inf - q)/tau
}

PROCEDURE rate(v (mV)) { :callable from hoc
        LOCAL a,q10
        q10=5^((celsius-23)/10)
        a = alp(v)
        inf = 1/(1 + a)
        tau = bet(v)/(q10*a0*(1+a))
}
