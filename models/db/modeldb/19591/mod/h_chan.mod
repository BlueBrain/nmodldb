TITLE h_chan.mod   squid sodium, potassium, and leak channels
 
COMMENT
%W%                            %G%
ENDCOMMENT
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX h_chan
        RANGE gbar, i
        NONSPECIFIC_CURRENT i
        GLOBAL ena, ek, minf, 
               am, bm, cm, dm, taum_min
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        dt (ms)
        gbar = 1.648e-4 (mho/cm2)
        ena = 55 (mV)
        ek = -90 (mV)

        am = 0.001
        bm = -0.1818
        cm = 0.5
        dm = -75
        taum_min = 1e-3

}
 
STATE {
        m
}
 
ASSIGNED {
        i (mA/cm2)
        minf
}
 
LOCAL mexp
 
BREAKPOINT {
        SOLVE states
        i = gbar*m*(v - ek)*2/3 +  gbar*m*(v - ena)/3
}
 
UNITSOFF
 
INITIAL {
	rates(v)
	m = minf
}

PROCEDURE states() {  :Computes state variables m, h, and n 
        rates(v)      :             at the current v and dt.
        m = m + mexp*(minf-m)
        VERBATIM
        return 0;
        ENDVERBATIM
}
 
PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
        LOCAL  tau,alpha,beta
        TABLE minf, mexp DEPEND dt FROM -100 TO 100 WITH 2000

                :"m" h current activation system
        alpha = am*exp(bm*cm*(v-dm))
        beta = am*exp(-bm*(1-cm)*(v-dm))
        tau = 1/(alpha + beta)
        minf = alpha*tau
        if (tau<taum_min) { tau = taum_min }
        mexp = 1 - exp(-dt/tau)

}
 
 
UNITSON

