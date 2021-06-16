TITLE BKDNaDR.mod   CNS sodium, potassium, calcium and leak channels
 
COMMENT
Active channels for a cortical pyramidal cell.
Fast Na and potassium (delayed rectifier); leak current
Specifications in Bernander, Koch and Douglas, J.Neurophys
72:2743-7253, 1994
BPG 2-4-97
ENDCOMMENT
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX BKDNaDR
        USEION na READ ena WRITE ina
        USEION k READ ek WRITE ik
        NONSPECIFIC_CURRENT il
        RANGE gnabar, gkbar, gl, el, ina, ik, il
        GLOBAL mss, hss, nss, mexp, hexp, nexp
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        dt (ms)
        gnabar = 0.2 (mho/cm2)
        ena = 50 (mV)
        gkbar = 0.12 (mho/cm2)
        ek = -95 (mV)
        gl = 0.00001 (mho/cm2)
        el = -66 (mV)
}
 
STATE {
        m h n
}
 
ASSIGNED {
        ina (mA/cm2)
        ik (mA/cm2)
        il (mA/cm2)
        mss hss nss mexp hexp nexp
}
 
BREAKPOINT {
        SOLVE states
        ina = gnabar*m*m*h*(v - ena)
        ik = gkbar*n*n*(v - ek)      
        il = gl*(v - el)
}
 
UNITSOFF
 
INITIAL {
	rates(v)
	m = mss
	h = hss
	n = nss
}

PROCEDURE states() {  :Computes state variables m, h, and n 
        rates(v)      :             at the current v and dt.
        m = m + mexp*(mss-m)
        h = h + hexp*(hss-h)
        n = n + nexp*(nss-n)
}
 
PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize ss at resting v.
        TABLE mss, mexp, hss, hexp, nss, nexp DEPEND dt FROM -100 TO 100 WITH 200
                :"m" sodium activation system
        mss = 1/(1+exp(-(v+40)/3))
        mexp = 1 - exp(-dt/0.05)
                :"h" sodium inactivation system
        hss = 1/(1+exp((v+45)/3))
        hexp = 1 - exp(-dt/0.5)
                :"n" potassium activation system
        nss = 1/(1+exp(-(v+40)/3))
        nexp = 1 - exp(-dt/2)
}
 
 
UNITSON

