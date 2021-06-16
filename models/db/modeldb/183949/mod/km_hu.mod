TITLE Borg-Graham like K-M channel, taken from Hu et al. J Neurosci 29:14472,2009

NEURON {
	SUFFIX km_hu
	USEION k READ ek WRITE ik
        RANGE gbar,inf,tau,m,ik
        GLOBAL vhalf,a0,zeta,gm
}
UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
        (mM) = (milli/liter)

}

PARAMETER {
	v (mV)
        ek (mV)
	tm=1.0
	celsius 	(degC)
	gbar=.0012 (mho/cm2)
        vhalf=-43   (mV)
        a0=0.004      (/ms)
        zeta=3.5  
        gm=0.5  
}



STATE {
        m
}

ASSIGNED {
	ik (mA/cm2)
        inf
        tau
}

INITIAL {
        rate(v)  
        m=inf
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	ik = gbar*m*(v-ek)

}

FUNCTION alp(v(mV)) {
  alp = a0*exp(1e-3*zeta*gm*(v-vhalf)*9.648e4/(8.315*(273.16+celsius)))
}

FUNCTION bet(v(mV)) {
  bet = a0*exp(-1e-3*zeta*(1-gm)*(v-vhalf)*9.648e4/(8.315*(273.16+celsius))) 
}

DERIVATIVE state {  
        rate(v)
        m' = (inf - m)/tau
}

PROCEDURE rate(v (mV)) { :callable from hoc
        LOCAL q10
        q10=5^((celsius-32)/10)
         inf = alp(v)/(alp(v)+bet(v))
        tau =(1/q10)*(tm+1/(alp(v)+bet(v)))
}














