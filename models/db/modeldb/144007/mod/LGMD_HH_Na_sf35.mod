TITLE HH_Na35 channel for LGMD
: Altered by PWJ, Jan 2011 to add malpha_sf to speed up the activation constant

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
}

NEURON {
    SUFFIX HH_Na35
    USEION na READ ena WRITE ina
    RANGE gmax
}

PARAMETER {
    gmax= 0.045 (mho/cm2)
}

ASSIGNED { 
    v (mV)
    ena (mV)
    ina (mA/cm2)
    malpha (/ms)
    mbeta (/ms)
    halpha (/ms)
    hbeta (/ms)
}

STATE {
    m h
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ina  = gmax*m*m*m*h*(v-ena)
}

INITIAL {
    settables(v)
    m = malpha/(malpha+mbeta)
    h = halpha/(halpha+hbeta)
}

DERIVATIVE states {  
    settables(v)      
    m' = 4*((malpha*(1-m)) - (mbeta*m))
    h' = 4*((halpha*(1-h)) - (hbeta*h))
}

UNITSOFF

PROCEDURE settables(v (mV)) {
    LOCAL den, den2, m_sf, h_sf, malpha_sf, mpos, hpos
    TABLE malpha, mbeta, halpha, hbeta
          FROM -100 TO 100 WITH 2000

    malpha_sf = 1
    m_sf = 3.5
    h_sf = 3.5
    mpos = 35
    hpos = 52

    den = (exp(-0.1*(v+mpos))-1) 
    if (v > (-mpos-.5) && v < (-mpos+.5)) {
      malpha = malpha_sf*m_sf*1
    } else {
      malpha = malpha_sf*m_sf*(-0.1*(v+mpos))/den
    }
    mbeta  = m_sf*4*exp(-(v+mpos+25)/12)

    halpha = h_sf*0.07*exp(-0.1*(v+hpos))
    den2 = (exp(-0.1*(v+hpos-30))+1)
    hbeta  = h_sf*1/den2
}

UNITSON


