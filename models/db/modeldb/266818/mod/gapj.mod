NEURON {
    POINT_PROCESS HalfGap
    POINTER vgap
    NONSPECIFIC_CURRENT i
    RANGE g, i, q10, vcell, vvgap
}

UNITS {
    (nA) = (nanoamp)
    (mV) = (millivolt)
    (pS) = (picosiemens)
}

PARAMETER {
    g = 1 (pS)
}

ASSIGNED {
    v (mV)
    celsius (degC)
    vgap (mV)
    i (nA)
    vvgap (mV)
    vcell (mV)
    q10 (1)
}

BREAKPOINT {
    vvgap = vgap
    vcell = v
    q10 = 1.43 ^((celsius-37(degC))/10 (degC))
    i = 1e-6 * g * q10 * (v - vgap)
}

:the 1e-6 is to fix the units to nA