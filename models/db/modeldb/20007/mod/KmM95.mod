:Migliore file Modify by Maciej Lazarewicz (mailto:mlazarew@gmu.edu) May/16/2001

TITLE Borg-Graham K-M channel

NEURON {
	SUFFIX KmM95
	USEION k READ ek WRITE ik
        RANGE  gbar,ik
        GLOBAL inf,tau
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
        (mM) = (milli/liter)
}

PARAMETER {
        cai 		(mM)
	v 		(mV)
        ek 	   	(mV)
	celsius = 30	(degC)
	gbar	= .003 	(mho/cm2)
        vhalf	= -55   (mV)
        a0	= 0.006 (/ms)
        zeta	= -10    (1)
        gm	= 0.06   (1)
        st	= 1
}

STATE { m }

ASSIGNED {
	ik (mA/cm2)
        inf
        tau
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	ik = gbar*m^st*(v-ek)
}

INITIAL {
	rate(v)
	m=inf
}

FUNCTION alp(v(mV)) {
  alp = exp( 1.e-3*zeta*(v-vhalf)*9.648e4/(8.315*(273.16+celsius)))
}

FUNCTION bet(v(mV)) {
  bet = exp(1.e-3*zeta*gm*(v-vhalf)*9.648e4/(8.315*(273.16+celsius))) 
}

DERIVATIVE state {
        rate(v)
        m' = (inf - m)/tau
}

PROCEDURE rate(v (mV)) { :callable from hoc
        LOCAL a,q10
        q10 = 5^((celsius-23)/10)

        a   = alp(v)
        inf = 1/(1 + a)
        tau = bet(v)/(q10*a0*(1+a))
}














