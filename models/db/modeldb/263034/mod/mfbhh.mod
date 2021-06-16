TITLE: mfbhh.mod    Sodium and potassium channels of mossy fiber boutons

COMMENT
  This is the Hodgkin-Huxley treatment for the set of sodium, potassium, 
  and leakage channels found in the hippocampal mossy fiber boutons.
  ("Presynaptic action potential amplification by voltage-gated Na+ channels in 
  hippocampal mossy fiber boutons" Neuron 45:405-417 (2005).)
  Global activation & inactivation shift; make vShift (Donnan) global by 12 mV.  
  "Engel & Jonas model (2005)" reconstructed by Kamiya 
ENDCOMMENT
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
	    (S) = (siemens)
}
 
? interface
NEURON {
        SUFFIX mfbhh
        USEION na READ ena WRITE ina
        USEION k READ ek WRITE ik
        NONSPECIFIC_CURRENT il
        RANGE gnabar, gkbar, gl, el, gna, gk
        GLOBAL minf, hinf, ninf, rinf, mtau, htau, ntau, rtau
	THREADSAFE : assigned GLOBALs will be per thread
}
 
PARAMETER {
        gnabar = 0.05 (S/cm2)
        gkbar = 0.036 (S/cm2)
        gl = .0001 (S/cm2)
        el = -81 (mV)
}
 
STATE {
        m h n r
}
 
ASSIGNED {
        v (mV)
        celsius (degC)
        ena (mV)
        ek (mV)

	    gna (S/cm2)
	    gk (S/cm2)
        ina (mA/cm2)
        ik (mA/cm2)
        il (mA/cm2)
        minf hinf ninf rinf
	mtau (ms) htau (ms) ntau (ms) rtau (ms)
}
 
? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
    gna = gnabar*m*m*m*h
	ina = gna*(v - ena)
    gk = gkbar*n*n*n*n*r
	ik = gk*(v - ek)      
    il = gl*(v - el)
}
 
 
INITIAL {
	rates(v)
	m = minf
	h = hinf
	n = ninf
    r = rinf
}

? states
DERIVATIVE states {  
        rates(v)
        m' = (minf-m)/mtau
        h' = (hinf-h)/htau
        n' = (ninf-n)/ntau
        r' = (rinf-r)/rtau
}
 
:LOCAL q10


? rates
PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
        LOCAL  alpha, beta, sum, q10
        TABLE minf, mtau, hinf, htau, ninf, ntau, rinf, rtau DEPEND celsius FROM -100 TO 100 WITH 200

UNITSOFF
        q10 = 3^((celsius - 23)/10)
     
	 :"m" sodium activation system
        alpha = 93.8285*vtrap(-(v-12-105.023),17.7094)
        beta =  0.168396*exp(-(v-12)/23.2707)
        sum = alpha + beta
	    mtau = 1/(q10*sum)
        minf = alpha/sum
     :"h" sodium inactivation system
        alpha = 0.000354*exp(-(v-12)/18.706)
        beta = 6.62694/(exp(-(v-12+17.6769)/13.3097)+1)
        sum = alpha + beta
	    htau = 1/(q10*sum)
        hinf = alpha/sum
		
     :"n" potassium activation system
        alpha = .01*vtrap(-(v+55),10) 
        beta = .125*exp(-(v+65)/80)
	    sum = alpha + beta
        ntau = 1/(q10*sum)
        ninf = alpha/sum
     :"r" potassium inactivation system
        alpha = 0.0000256077*exp(-v/45.4217)
        beta = 0.0330402/(exp(-(v+45.6599)/2.30235)+1)  :Recombinant Kv1.4
        sum = alpha + beta
	    rtau = 1/(q10*sum)
        rinf = alpha/sum
}
 
FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON
