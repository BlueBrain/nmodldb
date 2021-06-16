TITLE leakdend.mod   leak channels
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX leakdend
	NONSPECIFIC_CURRENT il
	RANGE gl, el  

}

PARAMETER {
        v (mV)
        celsius (degC)
        dt (ms)
	gl = 0.000125 (mho/cm2)
        el = -30 (mV)
}
 
STATE {
        c
}
 
ASSIGNED {
 
        il (mA/cm2)
}
 
BREAKPOINT {
        il = gl*(v - el)
}

UNITSON

 



