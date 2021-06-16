TITLE P-type calcium channel

COMMENT

P/Q-type VGCC. Starting from the formulation for P/Q-type VGCCs in Anwar et al., 2012, we
multiplied the activation curve by a sigmoid function to account for the fact that we did not observe
P/Q channel activation below -50 mV. We also reduced the activation time by 40% to reproduce the
observed Ca 2+ spiking rate at depolarised states.

Current Model Reference: Karima Ait Ouares , Luiza Filipis , Alexandra Tzilivaki , Panayiota Poirazi , Marco Canepari (2018) Two distinct sets of Ca 2+ and K + channels 
are activated at different membrane potential by the climbing fibre synaptic potential in Purkinje neuron dendrites. 

PubMed link: 

Contact: Filipis Luiza (luiza.filipis@univ-grenoble-alpes.fr)

ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    SUFFIX newCaP1
    USEION ca READ cai, cao WRITE ica
    RANGE pcabar, ica, gk, vhalfm, cvm, vshift,pp, t1, t2
}

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (nA) = (nanoamp)
    (pA) = (picoamp)
    (S)  = (siemens)
    (nS) = (nanosiemens)
    (pS) = (picosiemens)
    (um) = (micron)
    (molar) = (1/liter)
    (mM) = (millimolar)     
}

CONSTANT {
    q10 = 3
    F = 9.6485e4 (coulombs)
    R = 8.3145 (joule/kelvin)
}

PARAMETER {
    v (mV)
    celsius (degC)

    cai (mM)
    cao (mM)
    pp=1
       t1=0.4
    t2=40

    vhalfm = -29.458 (mV)
    cvm = 8.429(mV)
    vhalfh = -11.039 (mV)
    cvh = 16.098 (mV)
    vshift = 0 (mV)

    pcabar = 0.00049568 (cm/s)
}

ASSIGNED {
    qt
    ica (mA/cm2)
    minf
    taum (ms)
    corr 
    gk (coulombs/cm3)
    T (kelvin)
    E (volt)
    zeta

}

STATE { m h }

INITIAL {
    qt = q10^((celsius-23 (degC))/10 (degC))
    T = kelvinfkt( celsius )
    rates(v)
    m = minf
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    
    ica = (1e3) * pcabar * m * m * m * gk
}

DERIVATIVE states {
    rates(v)
    m' = (minf-m)/taum
}

FUNCTION ghk( v (mV), ci (mM), co (mM), z )  (coulombs/cm3) { 
    E = (1e-3) * v
      zeta = (z*F*E)/(R*T)  
    
    if ( fabs(1-exp(-zeta)) < 1e-6 ) {
        ghk = (1e-6) * (z*F) * (ci - co*exp(-zeta)) * (1 + zeta/2)
    } else {
        ghk = (1e-6) * (z*zeta*F) * (ci - co*exp(-zeta)) / (1-exp(-zeta))
    }
}

PROCEDURE rates( v (mV) ) {
    corr=1/(1+exp(-t1*(v+t2)))

    minf =( 1 / ( 1 + exp(-(v-vhalfm-vshift)/cvm) ))*corr

    taum = pp*taumfkt(v-vshift)/qt
    
    gk = ghk(v-vshift, cai, cao, 2)
}


FUNCTION kelvinfkt( t (degC) )  (kelvin) {
    UNITSOFF
    kelvinfkt = 273.19 + t
    UNITSON
}

FUNCTION taumfkt( v (mV) ) (ms) {
    UNITSOFF
    if (v>=-40) {
        taumfkt = 0.2702 + 1.1622 * exp(-(v+26.798)*(v+26.798)/164.19)
    } else {
        taumfkt = 0.6923 * exp(v/1089.372)
    }
    UNITSON
}

