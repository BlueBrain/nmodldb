TITLE Leak current
: FORREST MD (2014) Two Compartment Model of the Cerebellar Purkinje Neuron
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX leak
        NONSPECIFIC_CURRENT il
        RANGE  gl, el
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        celsius = 37 (degC)
        dt (ms)
        gl = .0003 (mho/cm2)
        el = -68 (mV)
}
  
ASSIGNED {
        il (mA/cm2)
}



BREAKPOINT {
 il = gl*(v - el)
}


