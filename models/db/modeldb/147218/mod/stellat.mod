: hyperpolarizing current due to activation of stellate cells

NEURON {
    SUFFIX stellat
    NONSPECIFIC_CURRENT  i
    RANGE i, e, g
}

PARAMETER {
    g = 2e-6 (siemens/cm2)  < 0, 1e9 >
    e = -95    (millivolts)
}

ASSIGNED {
    i   (milliamp/cm2)
    v   (millivolt) 
}

BREAKPOINT { i = g*(v - e) }
