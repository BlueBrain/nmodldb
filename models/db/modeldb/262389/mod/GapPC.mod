TITLE Generic gap junctions compatible with parallel context
: Gap junctions


NEURON {
    POINT_PROCESS GapPC
    RANGE r, i, vgap
    ELECTRODE_CURRENT i
}

PARAMETER {
    r = 1428 (megohm)  : 0.7 nS for AII, Veruki ML & Hartveit E (2008)
}

ASSIGNED {
    v     (mV)
    vgap  (mV)
    i     (nA)
}

BREAKPOINT {
    i = (vgap - v) / r
}

