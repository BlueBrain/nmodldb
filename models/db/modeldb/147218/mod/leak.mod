: leak current

NEURON {
    SUFFIX leak
    NONSPECIFIC_CURRENT  i
    RANGE i, e, g
}

PARAMETER {
    g = 2e-5 (siemens/cm2)  < 0, 1e9 >
    e = -60 (millivolts)
}

ASSIGNED {
    i   (milliamp/cm2)
    v   (millivolt) 
}

BREAKPOINT { i = g*(v - e) }
