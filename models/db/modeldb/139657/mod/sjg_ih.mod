TITLE sjg_ih.mod  Ih conductances

COMMENT

Ih for VCN neurons - average from several studies in auditory neurons
Edits by SJG and MHH
Based on the Leao MNTB model (Leao et al. 2005)

fit to MNTB data by Sarah J. Griffin, MRC Toxicology Unit, Leicester

ENDCOMMENT

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (nA) = (nanoamp)
}

NEURON {
    SUFFIX sjg_ih
    NONSPECIFIC_CURRENT i
    RANGE ghbar, gh, ih, sjgeh
    GLOBAL uinf, utau
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
    v (mV)
    celsius = 22 (degC)
    dt (ms)
    ghbar = 0.00318 (mho/cm2) <0,1e9>
    sjgeh = -43 (mV)
}

STATE {
    u
}

ASSIGNED {
    gh (mho/cm2)
    i (mA/cm2)
    uinf
    utau (ms)
}

LOCAL uexp

BREAKPOINT {
    SOLVE states
    
    gh = ghbar*u
    i = gh*(v - sjgeh)
}

UNITSOFF

INITIAL {
    trates(v)
    u = uinf
}

PROCEDURE states() {  :Computes state variables m, h, and n
    trates(v)      :             at the current v and dt.
    u = u + uexp*(uinf-u)
VERBATIM
    return 0;
ENDVERBATIM
}

LOCAL q10
PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
		      
    q10 = 3^((celsius - 22)/10)
    uinf = 1 / (1+exp((v + 101) / 11))
    utau = (10000 / (235.55*exp(0.0782*(v+23.76)) + 0.33*exp(-0.0614*(v+23.76)))) + 154.57

}

PROCEDURE trates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
    LOCAL tinc
    TABLE uinf, uexp
    DEPEND dt, celsius FROM -200 TO 150 WITH 350

    rates(v)    : not consistently executed from here if usetable_hh == 1
        : so don't expect the tau values to be tracking along with
        : the inf values in hoc

    tinc = -dt * q10
    uexp = 1 - exp(tinc/utau)
}

FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
    if (fabs(x/y) < 1e-6) {
        vtrap = y*(1 - x/y/2)
    }else{
        vtrap = x/(exp(x/y) - 1)
    }
}

UNITSON
