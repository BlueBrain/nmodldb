TITLE Cone photoreceptor Kx h Ca channel
: Cone photoreceptor Kx h Ca channel
:
: Based on Publio et al. (2009)


NEURON {
    SUFFIX ConePR
    USEION Ca WRITE iCa VALENCE 2
    USEION Cl WRITE iCl  VALENCE 1
    USEION Kca WRITE iKca VALENCE 1
    
    NONSPECIFIC_CURRENT il, iCGMP
    
    RANGE gCabar, gCa, eCa, SCa, VhalfCa, aoCa
    
    RANGE gClbar,gCl, eCl, SCl
    RANGE gKcabar,gKca, eKca
    
    RANGE gl, el
    RANGE gCGMP, eCGMP
    
    :temporal parameters
    RANGE FactorCaI
    RANGE mCl, Cas
}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (mS) = (millimho)
    (mol)= (1)
    (M)  = (mol/liter)
    (uM) = (micro M)
}

PARAMETER {
    : Calcium channel 
    gCabar  = 4.92 (mS/cm2) < 0, 1e9 >
    eCa     = 40 (mV)
    aoCa    = 0.0031 (/ms)
    VhalfCa = -16.6 (mV)
    SCa     = 5.7   (mV)

    : Cl channel    
    eCl       = -45 (mV)
    gClbar    = 6.5 (mS/cm2) < 0 , 1e9 >
    SCl       = 0.09 (uM)
    Clh       = 0.37 (uM)
    FactorCaI = 0.45
 
    : Ca-dependent K current
    eKca    = -80 (mV)
    gKcabar = 0.5 (mS/cm2)
 
    : leak
    gl = 0.01 (mS/cm2)
    el = 0 (mV)

    : cGMP gated channel 
    gCGMP= 0   (mS/cm2)
    :1.8   (mS/cm2)
    eCGMP=0.8 (mV)
}

STATE {
    nCa
    mKca
}

ASSIGNED {
    v     (mV)
    iCa   (mA/cm2)
    il    (mA/cm2)
    iCl   (mA/cm2)
    iCGMP (mA/cm2)
    iKca  (mA/cm2)

    : Ca-dependent potassium channel, Kca
    infmKca
    taumKca (ms)
    
    infCa
    tauCa  (ms) 
    
    Cas  (uM)
    mCl

    : the paremeter for activation
    mKca1
    gKca (mho/cm2)
    gCa (mho/cm2)
    gCl (mho/cm2)
}

INITIAL {
    rate(v)
    nCa  = infCa
    mKca = infmKca
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    gCa = (0.001)*gCabar*nCa
    iCa = gCa*(v - eCa)
    
    UNITSOFF
    :if (iCa >= 0) 
    :{
    :    Cas =0
    :}
    :if (iCa < 0) 
    :{
        Cas =-0.2+FactorCaI * (-iCa) * 1 *  0.5         /(1.6e-19)/  (6.023e23) * 1e-6         *1e14    
    :                  mA/cm2 * ms-> n coul/cm2  ->n e /cm2-> nmol/cm2  -> mol /cm2     scale factor
    : all the calculation without consideration of volume
         :    }
    
    mCl = 1/(1+ exp ( (Clh - Cas)/ SCl  ) ) 
    gCl = (0.001)* gClbar * mCl
    iCl = gCl*(v-eCl)   
    
    mKca1=Cas/(Cas+0.3)
    gKca=(0.001)*gKcabar*mKca*mKca*mKca1
    iKca=gKca*(v-eKca) 
    
    UNITSON
    
    il  = (0.001)*gl*(v-el)
      
    iCGMP = (0.001)*gCGMP*(v-eCGMP)
    
    : the current is in the unit of mA/cm2
}

DERIVATIVE states {
    rate(v)
    nCa' = (infCa - nCa)/tauCa
    mKca'= (infmKca - mKca ) /taumKca
}

UNITSOFF

FUNCTION alphamKca(v(mV)) (/ms) { 
    alphamKca = (0.001)*15*(80-v)/ ( exp( (80-v)/40 ) -1)
    :alter from orginal settings where it is in the unit of 1/s
}

FUNCTION  betamKca (v(mV)) (/ms) {
    betamKca = (0.001)*20*exp (-v/35)
}

UNITSON


FUNCTION alphaCa(v(mV))(/ms) { 
    alphaCa = aoCa*exp( (v - VhalfCa)/(2*SCa)   )
}

FUNCTION betaCa(v(mV))(/ms) { 
    betaCa = aoCa*exp( - ( v-VhalfCa)/(2*SCa) )
}

PROCEDURE rate(v (mV)) { LOCAL a, b
    a = alphamKca(v)
    b = betamKca(v)
    taumKca = 1/(a + b)
    infmKca = a/(a + b)
    
    a = alphaCa(v)
    b = betaCa(v)
    tauCa = 1/(a + b)
    infCa = a/(a + b)
}

