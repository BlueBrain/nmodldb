TITLE KCa channels and Ca mechanism for LGMD SFA

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
}

NEURON {
    SUFFIX KCa
    USEION ca READ cai
    USEION k READ ek WRITE ik
    RANGE gmax, kD_ca, ikca, cm
}

PARAMETER {
    gmax = 0.04 (moh/cm2)
    kD_ca = .03 (mM)
}

ASSIGNED {
    v (mV)
    cai (mM)
    ek (mV)
    ik (mA/cm2)
    ikca (mA/cm2)
    cm (1) 
}

BREAKPOINT {
    cm = cai/(cai+kD_ca)
    ikca = gmax*cm*(v-ek)
    ik = ikca
}
