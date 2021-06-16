TITLE Kdr_chan.mod   delayed rectifier, granule cell. 
 
COMMENT
%W%                                 %G%
ENDCOMMENT
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX kdr_chan
        RANGE gbar, i
        NONSPECIFIC_CURRENT i
        GLOBAL e, minf, hinf,
               am, bm, cm, dm, taum_min,
               a1h, b1h, c1h, d1h, e1h, a2h, b2h, c2h, tauh_min
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        dt (ms)
        e = -90 (mV)

        gbar = 1.939e-2 (mho/cm2)

        am = 0.17
        bm = 0.091
        cm = 0.8
        dm = -38
        taum_min = 1e-3

        a1h = 0.76e-3
        b1h = 0.7e-3
        c1h = 6e-5
        d1h = -0.08
        e1h = -46
        a2h = 110e-5
        b2h = -0.0807
        c2h = 78
        tauh_min = 150
}
 
STATE {
        m h
}
 
ASSIGNED {
        i (mA/cm2)
        minf hinf
}
 
LOCAL mexp, hexp
 
BREAKPOINT {
        SOLVE states
        i = gbar*m*m*m*m*h*(v - e)
}
 
UNITSOFF
 
INITIAL {
	rates(v)
	m = minf
	h = hinf
}

PROCEDURE states() {  :Computes state variables m, h, and n 
        rates(v)      :             at the current v and dt.
        m = m + mexp*(minf-m)
        h = h + hexp*(hinf-h)
        VERBATIM
        return 0;
        ENDVERBATIM
}
 
PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
        LOCAL  tau,alpha,beta
        TABLE minf, mexp, hinf, hexp DEPEND dt FROM -100 TO 100 WITH 2000

                :"m" sodium activation system
        alpha = am*exp(bm*cm*(v-dm))
        beta = am*exp(-bm*(1-cm)*(v-dm))
        tau = 1/(alpha + beta)
        minf = alpha*tau
        if (tau<taum_min) { tau = taum_min }
        mexp = 1 - exp(-dt/tau)

                :"h" sodium inactivation system
        if (v>e1h) {
          alpha = a1h
        }
        else {
          alpha = b1h+c1h*exp(d1h*(v-e1h))
        }
        beta = a2h/(1+exp(b2h*(v+c2h)))
        tau = 1/(alpha + beta)
        hinf = alpha*tau
        if (tau<tauh_min) { tau = tauh_min }
        hexp = 1 - exp(-dt/tau)
}
 
 
UNITSON

