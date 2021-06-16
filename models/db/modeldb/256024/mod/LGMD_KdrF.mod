TITLE KdrF channel for LGMD
: RBD

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
	(S) = (siemens)
}

NEURON {
    THREADSAFE
    SUFFIX KdrF
    USEION k READ ek WRITE ik
    RANGE gmax, g, t1, t2	: Range variables can differ across neuron (non constants)
    GLOBAL vhalf
}

PARAMETER {
    gmax= 0.008 (mho/cm2)
    
    vhalf=-39	(mV)
	vn2=-40	(mV)
	vl=-55	(mV)
	t1=2.7	(ms)
	t2=0.1	(ms)
	tns=7	(mV)
	zn=9	(mV)
	zl=-9	(mV)
}

ASSIGNED { 
    v (mV)
    ek (mV)
    
    ik (mA/cm2)
	ninf
    ntau (ms)
    g (S/cm2)
}

STATE {
    n
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    g  = gmax*n
    ik  = g*(v-ek)

}

INITIAL {
    settables(v)
    n = ninf
}

DERIVATIVE states {  
    settables(v)    
    n' = (ninf - n)/ntau
}

UNITSOFF

PROCEDURE settables(v (mV)) {
	TABLE ninf, ntau DEPEND vhalf, t1
          FROM -100 TO 50 WITH 600

	ninf = 1/(1 + exp((vhalf-v)/zn))^2
	ntau = 2*t1/(1+exp((vn2-v)/tns))/(1+ exp((vl-v)/zl))+t2
	:ntau = 4*t1/(1+exp((vn2-v)/tns))*ninf+t2

}

UNITSON


