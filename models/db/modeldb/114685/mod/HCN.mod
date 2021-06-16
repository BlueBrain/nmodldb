TITLE Potassium Ih channel for GPi neuron model

COMMENT

 HCN1 and HCN2 channels are present in the EP (see Fig4 of Chan2004)
 HCN2 immunoreactivity appeared to be a bit stronger, but this may be
 due to the methods used in the study.  Ih modeled from an LGN relay neuron:
    -Huguenard & McCormick 1992 (recordings at 35.5degC -- disregard Q10)
    -McCormick & Pape 1990 (recordings at 36.1degC -- disregard Q10)
 Code amended from Gillies2006

ENDCOMMENT

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    SUFFIX Ih
    NONSPECIFIC_CURRENT ih
    RANGE gh
    GLOBAL eih
}

PARAMETER {
    v (mV)
    dt (ms)
    gh    = 0.001 (mho/cm2)
    eih   = -56.11 (mV)
    celsius
}

STATE {
    f
}

ASSIGNED { 
    ih (mA/cm2)
    finf  
    ftau (ms)
}

BREAKPOINT {
    SOLVE integrate METHOD cnexp
    ih = gh*f*(v-eih)
}

UNITSOFF

INITIAL {
    setinf(v)
    f = finf
}

DERIVATIVE integrate {
    setinf(v)
    f' = (finf - f)/ftau
}

PROCEDURE setinf(v) {
    
    LOCAL vhalf, sfactor
    TABLE finf, ftau DEPEND celsius FROM -100 TO 100 WITH 400

    vhalf = -75     : membrane potential when Gh is half-activated (mV)
    sfactor = 5.5   : slope factor which determines the steepness of the fitted curve
    finf = 1.0/(1+exp((v-vhalf)/sfactor))
    ftau = (1.0/(exp(-14.59-0.086*v)+exp(-1.87+0.0701*v)))
}

UNITSON
