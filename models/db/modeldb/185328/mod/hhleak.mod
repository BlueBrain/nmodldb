
TITLE HH Leak channel

UNITS {
        (mV) = (millivolt)
        (mA) = (milliamp)
        (S) = (siemens)
}

NEURON {
        SUFFIX HHleak
        NONSPECIFIC_CURRENT i
        RANGE g, e
}

PARAMETER {
        g = .0003       (S/cm2) <0,1e9>
        e = -54.3       (mV)
}

ASSIGNED {v (mV)  i (mA/cm2)}

BREAKPOINT {
        i = g*(v - e)
}