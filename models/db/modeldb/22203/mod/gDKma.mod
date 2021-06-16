COMMENT

dummy K conductance

ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    SUFFIX gDKkin3
    USEION k READ ek WRITE ik
    RANGE  gk, gbar
    GLOBAL VoltageOffset, gv, tau, tstart
}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (pS) = (picosiemens)
    (um) = (micron)
}

PARAMETER {
        gv[300]                             (ps/um2)
        tau[300]                            (ms)        : the tau[V] for the phaen Exp Fit
        gbar            = 5                 (pS/um2)    : 0.03 mho/cm2
        v                                   (mV)
        VoltageOffset   = 100.5             (mV)        : but modified by fit.hoc
        tstart          = 10                (ms)        : later modified by fit.hoc
        dt                                  (ms)
        measTime                            (ms)        : Set by hoc to calculate factor
         celsius                             (degC)
}


ASSIGNED {
        ik      (mA/cm2)
        gk      (pS/um2)
        ek      (mV)
       }



BREAKPOINT {
        SOLVE states
        ik = (1e-4) * gk * (v - ek)
}


PROCEDURE states() {        
                          
        gbar = gv[v+VoltageOffset]
        if(t<tstart) { gk=0 }
        if (t>tstart) {gk = gbar}
        VERBATIM
        return 0;
        ENDVERBATIM
}



