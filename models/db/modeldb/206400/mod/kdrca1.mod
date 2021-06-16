TITLE K-DR channel
: from Klee Ficker and Heinemann
: modified to account for Dax et al.
: M.Migliore 1997
: MM: ntaufac added to model faster channel activation


UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)

}

PARAMETER {
    v (mV)
        ek (mV)     : must be explicitely def. in hoc
    celsius     (degC)
    gbar=.003 (mho/cm2)
   : ntaufac = 1 (1)      : factor that is multiplied with the expression for ntau

        vhalfn=13   (mV)    
        a0n=0.02      (/ms)
        zetan=-3    (1)
        gmn=0.7  (1)
    nmax=2  (1)
    q10=1
}

NEURON {
    SUFFIX kdr
    USEION k READ ek WRITE ik
        RANGE gkdr,gbar, ikdr, vhalfn
    GLOBAL ninf,taun
}

STATE {
    n
}

ASSIGNED {
    ik (mA/cm2)
    ikdr (mA/cm2)
        ninf
        gkdr
        :gbar
        taun
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    gkdr = gbar*n
    ik = gkdr*(v-ek)
    ikdr = gkdr*(v-ek)
    ik = ikdr

}

INITIAL {
    rates(v)    
    n=ninf
}


FUNCTION alpn(v) {
  alpn = exp(1.e-3*zetan*(v-vhalfn)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION betn(v) {
  betn = exp(1.e-3*zetan*gmn*(v-vhalfn)*9.648e4/(8.315*(273.16+celsius))) 
}

DERIVATIVE states {     : exact when v held constant; integrates over dt step
        rates(v) : ,vhalfn)
        n' = (ninf - n)/taun
}

PROCEDURE rates(v) { :callable from hoc
        LOCAL a,qt
        qt=q10^((celsius-24)/10)
        a = alpn(v)
        ninf = 1/(1+a)
        taun = betn(v)/(qt*a0n*(1+a))
    if (taun<nmax) {taun=nmax}
}














