TITLE hh.mod   squid sodium, potassium, and leak channels
 
COMMENT
 This is the original Hodgkin-Huxley treatment for the set of sodium, 
  potassium, and leakage channels found in the squid giant axon membrane.
The rate function parameters have been changed to correspond to
Dodge & Cooley (1973) "Action Potential of the Motorneuron"
IBM J. Res. Develop. May 219--229
ENDCOMMENT
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
	(S) = (siemens)
}
 
? interface
NEURON {
        SUFFIX dc
        USEION na READ ena WRITE ina
        USEION k READ ek WRITE ik
        NONSPECIFIC_CURRENT il
        RANGE gnabar, gkbar, gl, el, gna, gk, shift
        GLOBAL minf, hinf, ninf, mtau, htau, ntau, vrest
}
 
PARAMETER {
	: for node
        gnabar = .6 (S/cm2)	<0,1e9>
        gkbar = .1 (S/cm2)	<0,1e9>
        gl = .003 (S/cm2)	<0,1e9>
        el = -54.3 (mV)
	vrest = 0 (mV)
	shift = 0
}
 
STATE {
        m h n
}
 
ASSIGNED {
        v (mV)
        ena (mV)
        ek (mV)

	gna (S/cm2)
	gk (S/cm2)
        ina (mA/cm2)
        ik (mA/cm2)
        il (mA/cm2)
        minf hinf ninf
	mtau (ms) htau (ms) ntau (ms)
}
 
? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gna = gnabar*m*m*m*h
	ina = gna*(v - ena)
        gk = gkbar*n*n*n*n
	ik = gk*(v - ek)      
        il = gl*(v - el)
}
 
 
INITIAL {
	rates(v - vrest - shift)
	m = minf
	h = hinf
	n = ninf
}

? states
DERIVATIVE states {  
        rates(v - vrest - shift)
        m' =  (minf-m)/mtau
        h' = (hinf-h)/htau
        n' = (ninf-n)/ntau
}
 
? rates
PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
        LOCAL  alpha, beta, sum

UNITSOFF
                :"m" sodium activation system
        alpha = .4 * vtrap(25 - v, 5)
:        beta =  .4 * vtrap(v - 55, 5)
:	 typo in paper; personal communication from F. Dodge.
        beta =  .4 * vtrap(v - 45, 5)
        sum = alpha + beta
	mtau = 1/sum
        minf = alpha/sum
                :"h" sodium inactivation system
        alpha = .28 * exp((10 - v)/20)
        beta = 4 / (exp((40 - v)/10) + 1)
        sum = alpha + beta
	htau = 1/sum
        hinf = alpha/sum
                :"n" potassium activation system
:        alpha = .2*vtrap(20 - v,10) 
:	 typo in paper; personal communication from F. Dodge.
        alpha = .02*vtrap(20 - v,10) 
        beta = .25*exp((10 - v)/80)
	sum = alpha + beta
        ntau = 1/sum
        ninf = alpha/sum
}
 
FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON
