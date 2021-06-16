TITLE TRPM8 current plus calcium adaptation
: based on Voets 2004 with some modifications to resemble data from Malkia 2007
: Written by Orio, P. & Olivares E.  - December 2014
:

NEURON {
    SUFFIX trpm8
    NONSPECIFIC_CURRENT im8
    RANGE caM8, DVinf, vhalf, DV
    RANGE p_ca, accel
    RANGE em8, gm8, am8, C, z
}


UNITS {
    R = (k-mole) (joule/degC)
    (mA) = (milliamp)
    (mV) = (millivolt)
    (mol) = (1)
    (molar) = (1/liter)
    (mM) = (millimolar)
} 

CONSTANT {
    F = 96500        (coulomb)        : moles do not appear in units
}


PARAMETER {
    gm8 =0.0035        (mho/cm2)
    dE = 9e3        (joule)
    C = 67
    z = 0.65

    em8 = 0            (mV)
    DVmin = 0        (mV)
    DVmax = 200            (mV)
    Kca = 0.0005        (mM)

    p_ca=0.01    
    tauca=15000        (ms)
    taudv = 80000    (ms)
    d = 1           (micron)
    accel = 1
    n = 1
}


STATE {
    caM8        (mM)
    DV            (mV)
}

INITIAL {
    caM8=0
    DV= DVinf
}

ASSIGNED {
    celsius    (degC)
    v       (mV)
    im8        (mA/cm2)
    vhalf    (mV)
    am8
    DVinf        (mV)
    
}

BREAKPOINT {
    SOLVE states METHOD cnexp

    vhalf=(1000)*(C*R*celsius - dE)/(z*F)+DV
    
    am8=1/(1+exp(-z*F*(v-vhalf)/((1000)*R*(celsius+273.15))))
    im8=gm8*am8*(v-em8)
}

DERIVATIVE states {
    caM8' = accel*( - p_ca * (10000) * im8 / (2 * F * d) - caM8 / tauca)
    DVinf = DVmin+(DVmax-DVmin)*(caM8)/(Kca+caM8)
    DV' = accel*(DVinf-DV)/taudv
}
