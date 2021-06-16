: depolarizing current due to tonic activation of parallel fibers

NEURON {
    SUFFIX pf
    NONSPECIFIC_CURRENT  i
    RANGE i, e, g
}

PARAMETER {
    g = 0 (siemens/cm2)  < 0, 1e9 >
    e = 0    (millivolts)
}

ASSIGNED {
    i   (milliamp/cm2)
    v   (millivolt) 
}

BREAKPOINT { i = g*(v - e) }
