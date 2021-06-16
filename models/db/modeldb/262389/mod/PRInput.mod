TITLE Photoreceptor input current for horizontal cells
: Simulated photoreceptor EPSC for horizontal cells


NEURON {
    POINT_PROCESS PRInput
    RANGE delSteady, durSteady, ampSteady, i
    ELECTRODE_CURRENT i
}

UNITS {
    (nA) = (nanoamp)
}

PARAMETER {
    delSteady (ms)
    durSteady (ms)  < 0, 1e9 >
    ampSteady (nA)
}

ASSIGNED {
    i (nA)
}

BREAKPOINT {
    : purturbation at the following time points
    at_time(delSteady)
    at_time(delSteady + durSteady)

    if (t > delSteady && t < delSteady + durSteady) {
        i = ampSteady
    } else {
        i = 0
    }
}

