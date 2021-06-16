COMMENT

	Persistent sodium current from Av-Ron and Vidal 1999
	Implemented by C. Weaver, 2003

ENDCOMMENT

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {
	celsius 	(degC)
	gbar=.00005 (mho/cm2)
	vhp=-56 (mV)
	ap=0.075	(/mV)
	ptau=5 (ms)
        v       (mV)
        ena      (mV)
}


NEURON {
	SUFFIX nap
	USEION na READ ena WRITE ina
        RANGE gbar,gnap
        GLOBAL inf,ptau
	RANGE tot	
}

STATE {
	p
}

ASSIGNED {
	ina (mA/cm2)
	tot (mA/cm2)
        gnap  (mho/cm2)
	inf
}

INITIAL {
        rate(v)
        p = inf
: printf( "nap ina=%g\n",ina)
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	gnap = gbar*p
	tot = gnap*(v - ena)
	ina = gnap*(v-ena)

}

FUNCTION expn(v (mV),a(/mV), vhalf(mV)) {
  	expn = exp(-2*a*(v-vhalf))
}

DERIVATIVE state {     : exact when v held constant; integrates over dt step
        rate(v)
        p' = (inf - p)/ptau
}

PROCEDURE rate(v (mV)) { :callable from hoc
	inf = 1/(1+expn(v,ap,vhp))
}















