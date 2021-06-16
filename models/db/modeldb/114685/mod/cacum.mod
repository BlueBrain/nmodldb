TITLE calcium accumulation for GPi neuron model

COMMENT 

    From: Gillies and Willshaw (2006) and as implemented by Miocinovic (2006)

    Calcium accumulation into a volume of area*depth next to the
    membrane with an exponential decay (time constant tau) to resting
    level (given by the global calcium variable cai0_ca_ion).

ENDCOMMENT

NEURON {
    SUFFIX cacum
    USEION ca READ ica WRITE cai
    GLOBAL depth, tau, cai0
}

UNITS {
    (mM) = (milli/liter)
    (mA) = (milliamp)
    F = (faraday) (coulombs)	: Faraday constant 
}

PARAMETER {
    depth = 200 (nm)      : assume volume = area * depth
    tau = 10 (ms)
    cai0 = 5e-5 (mM)
    cai0_ca_ion
}

ASSIGNED {
    ica (mA/cm2)
}

STATE {
    cai (mM)
}

INITIAL {
    cai=cai0
}

BREAKPOINT {
    SOLVE integrate METHOD derivimplicit
}

DERIVATIVE integrate {
    cai' = -ica/depth/F/2 * (1e7) + (cai0 - cai)/tau
}