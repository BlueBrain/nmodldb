TITLE Na_chan.mod  na channel, granule cell. 
 
COMMENT
%W%                                 %G%
ENDCOMMENT
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX naf_chan
        RANGE gbar, i
        NONSPECIFIC_CURRENT i
        GLOBAL e, minf, hinf,
               am, bm, cm, dm, taum_min,
               ah, bh, ch, dh, tauh_min
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        dt (ms)
        e = 55 (mV)

        gbar = 7.43e-2 (mho/cm2)

        am = 1.5
        bm = 0.147
        cm = 0.55
        dm = -39
        taum_min = 0.05

        ah = 0.12
        bh = -0.177548
        ch = 0.5
        dh = -50
        tauh_min = 0.225
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
        i = gbar*m*m*m*h*(v - e)
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
        alpha = ah*exp(bh*ch*(v-dh))
        beta = ah*exp(-bh*(1-ch)*(v-dh))
        tau = 1/(alpha + beta)
        hinf = alpha*tau
        if (tau<tauh_min) { tau = tauh_min }
        hexp = 1 - exp(-dt/tau)
}
 
 
UNITSON

