: Calcium ion accumulation with radial and longitudinal diffusion
NEURON {
    SUFFIX jcal
    USEION ca READ cai, ica WRITE cai:, ica
    RANGE vrat
}

UNITS {
    (molar) = (1/liter)
    (mM) = (millimolar)
    (um) = (micron)
    (mA) = (milliamp)
    FARADAY = (faraday) (10000 coulomb)
    PI = (pi) (1)
}

PARAMETER {
    k1buf = 3 (/mM-ms) : Yamada et al. 1989
    k2buf = 0.005 (/ms)
    TotalBuffer = 5e-3 (mM)
    depth = 1 (um)
    taur = 1e1 (ms)	
    cainf = 5e-5 (mM)
    cain = 0.5 (1)
    gcaca = 0.2e-4
    cacamid = 1e-4
    cacas = 2e-5
}

ASSIGNED {
    ica (mM)
    cai (mM)
    vrat (1)
    Kd (/mM)
    B0 (mM)
}

STATE {
    ca (mM)
    cabuffer (mM)
    buffer (mM)
}

BREAKPOINT { 
    SOLVE state METHOD sparse
}

INITIAL {
    cai = cainf
    ca = cainf
    Kd = k1buf/k2buf
    B0 = TotalBuffer/(1 + Kd*cai)
    buffer = B0
    cabuffer = TotalBuffer-B0
    vrat = PI*(0.5-(0.5/(10.9495-1)/2)/2)*2*(0.5/(10.9495-1)/2)
}

LOCAL var1, var2, var3, dsqvol

KINETIC state {
    COMPARTMENT diam*diam*vrat {ca cabuffer buffer}
    var1 = (cainf-cai)/(0.35*taur)
    var2 = -cain*ica/(FARADAY*depth)
    var3 = gcaca/(1+exp((cacamid-ca)/cacas))
    ~ ca << (var1 + var2 + var3)
    dsqvol = diam*diam*vrat
    ~ ca + buffer <-> cabuffer (k1buf*dsqvol, k2buf*dsqvol)
    cai = ca
}
