TITLE Borg-Graham nonsecific cation H channel

COMMENT
Extended Hodgkin-Huxley model.
From Lyle J. Borg-Graham, Interpretations of data and mechanisms for
hippocampal pyramidal cell models.  In "Cerebral Cortex, Vol 13:
Cortical Models", Plenum Press 1998

Implemented by BPG 24-2-99
Renamed and tidied up BPG 3-11-99
ENDCOMMENT

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)

}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    SUFFIX L98NCh
    NONSPECIFIC_CURRENT ih
        RANGE ghbar,gh,ih
        GLOBAL ninf,linf,taul,taun,alphan,betan,alphal,betal,t0l
}

PARAMETER {
        dt (ms)
    v (mV)
        eh=-17 (mV)
    celsius = 32    (degC)
    ghbar=7.8e-5 (mho/cm2)
        vhalfn=0   (mV)
        vhalfl=-98   (mV)
        zn=0    (1)
        zl=-2    (1)
        gmn=0   (1)
        gml=0   (1)
    Kn=0    (1)
    Kl=0    (1)
    t0n=0   (ms)
    t0l=180 (ms)
}

STATE {
        l
}

ASSIGNED {
    ih (mA/cm2)
        gh (mho/cm2)
        ninf
        linf
    facn
    facl      
        taul
        taun
    alphan
    betan
    alphal
    betal
}

INITIAL {
    rates(v)
    l = linf
}

BREAKPOINT {
    SOLVE states
    gh = ghbar*l
    ih = gh*(v-eh)
}


FUNCTION alphap(v(mV),K,z,gamma,vhalf) {
  alphap = K*exp(z*gamma*(v-vhalf)*96.487/(8.314*(273.16+32))) 
}

FUNCTION betap(v(mV),K,z,gamma,vhalf) {
  betap = K*exp(-z*(1-gamma)*(v-vhalf)*96.487/(8.314*(273.16+32))) 
}

PROCEDURE states() {     : exact when v held constant; integrates over dt step
        rates(v)
        l = l + facl*(linf - l)
}

PROCEDURE rates(v (mV)) { :callable from hoc
    LOCAL q10
        TABLE linf, facl, taul DEPEND dt, t0l, celsius FROM -100 TO 100 WITH 200

: inactivation
        linf = 1/(1+betap(v,1.0,zl,0.0,vhalfl))
    taul = t0l
        facl = (1 - exp(-dt/taul))
}
