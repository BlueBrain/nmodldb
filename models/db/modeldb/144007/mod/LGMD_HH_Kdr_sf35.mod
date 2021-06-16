TITLE HH_Kdr35 channel for LGMD

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
}

NEURON {
    SUFFIX HH_Kdr35
    USEION k READ ek WRITE ik
    RANGE gmax, gk
}

PARAMETER {
    gmax= 0.018 (mho/cm2)
}

ASSIGNED { 
    v (mV)
    ek (mV)
    ik (mA/cm2)
    nalpha (/ms)
    nbeta (/ms)
    gk (mS/cm2)
}

STATE {
    n
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    gk  = gmax*n*n*n*n
    ik  = gmax*n*n*n*n*(v-ek)
}

INITIAL {
    settables(v)
    n = nalpha/(nalpha+nbeta)
}

DERIVATIVE states {  
    settables(v)      
    n' = 4*((nalpha*(1-n)) - (nbeta*n))
}

UNITSOFF

PROCEDURE settables(v (mV)) {
    LOCAL den, n_sf, npos
    TABLE nalpha, nbeta
          FROM -100 TO 100 WITH 2000

    n_sf = 3.5
    npos = 34

    den = (exp(-0.1*(v+npos))-1)
    if (v > (-npos-.5) && v < (-npos+.5)) { 
      nalpha = n_sf*0.1
    } else {
      nalpha = n_sf*(-0.01*(v+npos))/den
    }
    nbeta  = n_sf*0.125*exp(-(v+npos+10)/25)
}

UNITSON


