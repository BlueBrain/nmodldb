TITLE passive membrane properties for GPi

COMMENT

Nakanishi1990 measured the average input resistances of EP neurons
  Type I  --> 126 MOhm (8.9286e-9 mho)
  Type II --> 145 MOhm (6.8966e-9 mho)
All recordings were at 36degC, so there is no need to add a temperature
correction factor.  Note that if simulation involves Cs+ or TEA, it is
necessary to adjust gpas accordingly (See Fig1D).

ENDCOMMENT

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    SUFFIX GPi
    RANGE gpas, epas
    NONSPECIFIC_CURRENT ipas
}

PARAMETER {
    v (mV)
    dt (ms)
    gpas  = 6.8027e-5 (mho/cm2) <0,1e9>     
    epas  = -60 (mV)                        
    celsius
}

ASSIGNED { 
    ipas (mA/cm2)
}

BREAKPOINT {
    ipas = gpas*(v - epas)
}
