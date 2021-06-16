TITLE Synapse for Cone and horizontal cells
: Synapse for cones and horizontal cells. Can be excitatory or inhibitory by
: setting the sign of the maximum conductance.


NEURON {
    POINT_PROCESS Synapse
    RANGE V_pre, v_th, v_slope, g_max, i
    NONSPECIFIC_CURRENT i
}

PARAMETER {
    v_th    = -30.88  (mV)
    v_slope = 10      (mV)
    g_max   = 0.250   (umho)
    : g_max   = 0.00256 (umho)  : maximal conductance dAMPA
}

ASSIGNED {
    V_pre  (mV)
    i      (nA)
    g      (nA)
}

BREAKPOINT {
    g = tanh( (V_pre - v_th) / v_slope ) + 1
    i = -g_max * 0.5 * g
}

