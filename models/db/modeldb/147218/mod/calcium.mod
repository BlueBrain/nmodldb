: Inner Ca calculation with longitudinal diffusion

NEURON {
     SUFFIX calcium
     USEION ca READ ica WRITE cai
     GLOBAL BT, Kd, k, cab, delt
}

UNITS {
     (molar)  = (1/liter)
     (mM)     = (millimolar)
     (um)     = (micron)
     (mA)     = (milliamp)
     FARADAY  = (faraday)    (10000 coulomb)
     PI       = (pi)         (1)
}

PARAMETER {
     Kd    = 0.001    (mM)      : buffer dissociation constant
     BT    = 0.15     (mM)      : total buffer concentration
     k     = 0.00001  (cm/ms)   : constant of pumping in the ER
     cab   = 0.00005  (mM)      : calcium concentration in inner dendritic core
     delt  = 0.3      (um)      : thickness of cytoplasmic compartment
}

ASSIGNED {
     diam     (um)
     ica      (mA/cm2)
}

STATE { cai (mM) }

BREAKPOINT { SOLVE state METHOD sparse }

LOCAL x, y

INITIAL {
    cai = 96e-6
}

KINETIC state {   
     COMPARTMENT PI*delt*(diam-delt)  {cai}
     ~ cai << ( (-(1/(1+(BT/Kd)/(1+cai/Kd)^2))*((ica*(diam/2))/(FARADAY*delt*(diam/2-delt))
         +(1e4)* ((2*k*(cai-cab)*(diam/2-delt))/(delt*(diam/2-delt))))*PI*delt*(diam/2-delt)) )
}




