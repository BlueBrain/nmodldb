TITLE Cone photoreceptor IKv
: Cone photoreceptor IKv
: 
: Based on Publio et al. (2009)


NEURON {
    SUFFIX ConeIKv
    USEION Kv WRITE iKv VALENCE 1
    RANGE gKv,gKvbar, eKv
}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (mS) = (millimho)
}

PARAMETER {
    : potassium rectifier
    gKvbar = 2.0 (mS/cm2) <0,1e9>
    eKv = -80 (mV)
}

STATE {
    mKv
    hKv
}

ASSIGNED {
    v (mV)
    iKv (mA/cm2)
              
    : potussium rectifier, K
    infmKv
    taumKv  (ms)
    infhKv
    tauhKv  (ms)
     
    gKv (mho/cm2)
}

INITIAL {
    rate(v)
    mKv = infmKv
    hKv = infhKv
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    gKv = (0.001)*gKvbar*mKv*mKv*mKv *hKv
    iKv = gKv*(v - eKv)
}

DERIVATIVE states {
    rate(v)
    mKv' = (infmKv - mKv)/taumKv
    hKv' =  (infhKv - hKv )/tauhKv
}

UNITSOFF

FUNCTION alphamKv(v(mV)) (/ms) {
    alphamKv = (0.001)*5*(100-v)/( exp( (100-v)/42) -1 )
    :alter from orginal settings where it is in the unit of 1/s
}

FUNCTION betamKv (v(mV)) (/ms) {
    betamKv = (0.001)*9*exp (- (v-20) /40 )
}

FUNCTION alphahKv (v(mV)) (/ms) {
    alphahKv = (0.001)*0.15 *exp (-v/22)
}

FUNCTION betahKv (v(mV)) (/ms) { 
    betahKv = (0.001)*0.4125/( exp ( ( 10-v)/7 ) +1 )
}

UNITSON

PROCEDURE rate(v (mV)) { LOCAL a, b

    a = alphamKv(v)
    b = betamKv(v)
    taumKv = 1/(a + b)
    infmKv = a/(a + b)
    
    a = alphahKv(v)
    b = betahKv(v)
    tauhKv = 1/(a + b)
    infhKv = a/(a + b)
}

