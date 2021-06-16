TITLE M99Ka.mod  - K (A)
 
COMMENT
From Migliore et al J. Comput. Neurosci. 7(1):5-15, 1999
Based on Hoffman et al (1997) CA1 data.
Activation time constant, forward and backward rates parameterized
to allow for differences between proximal and distal Ka.
BPG 31-10-99
Ability to shift (in)activation curves and change time constants
added via parameters.
BPG 12-11-99
Minimum inactivation time constant added
BPG 9-1-01
ENDCOMMENT
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX M99Ka
        USEION k WRITE ik
        RANGE gkbar,gk,ik,tmfac,afac,bfac,thfac,vms,vhs,minf,hinf,mexp,hexp,thmin
:        GLOBAL minf, hinf, mexp, hexp
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        celsius = 35 (degC)
        dt (ms)
        gkbar = 0.048 (mho/cm2)
        ek = -90 (mV)
        tmfac = 4 (1)   : activation time constant factor (def. proximal)
        thfac = 1 (1)   : inactivation time constant factor
        thmin = 2 (ms)  : minimum inactivation time constant
        afac = 1.5 (1)  : activation forward rate factor
        bfac = 0.825 (1)    : activation backward rate factor
    vms = 0 (mV)    : activation curve voltage shift (+ve=right)
    vhs = 0 (mV)    : inactivation curve voltage shift
}
 
STATE {
        m h
}
 
ASSIGNED {
        ik (mA/cm2)
        minf hinf mexp hexp
}
 
BREAKPOINT {
        SOLVE states
        ik = gkbar*m*h*(v - ek)
}
 
UNITSOFF
 
INITIAL {
    rates(v)
    m = minf
    h = hinf
}

PROCEDURE states() {  :Computes state variables m, h
        rates(v)      :             at the current v and dt.
        m = m + mexp*(minf-m)
        h = h + hexp*(hinf-h)
}
 
PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
        LOCAL  alpha, beta, sum, tau
        TABLE minf,mexp,hinf,hexp DEPEND dt,vms,vhs,tmfac,afac,bfac,thfac,thmin FROM -100 TO 100 WITH 200
    :"m" potassium activation system
        alpha = exp(-0.038*(afac+1/(1+exp((v-vms+40)/5)))*(v-vms-11))
        beta = exp(-0.038*(bfac+1/(1+exp((v-vms+40)/5)))*(v-vms-11))
        sum = 1 + alpha
        minf = 1/sum
        tau = tmfac*beta/sum
        if (tau < 0.1) {tau=0.1}
        mexp = 1 - exp(-dt/tau)
    :"h" potassium inactivation system
        alpha = exp(0.11*(v-vhs+56))
        hinf = 1/(1+alpha)
        tau = thfac*0.26*(v+50)
        if (tau < thmin) {tau=thmin}
        hexp = 1 - exp(-dt/tau)
}
 
UNITSON
