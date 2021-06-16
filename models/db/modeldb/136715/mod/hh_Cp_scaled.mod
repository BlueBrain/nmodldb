TITLE hh_Cp_scaled.mod   squid sodium, potassium, and leak channels
 
COMMENT
 This is the original Hodgkin-Huxley treatment for the set of sodium, 
  potassium, and leakage channels found in the squid giant axon membrane.
  ("A quantitative description of membrane current and its application 
  conduction and excitation in nerve" J.Physiol. (Lond.) 117:500-544 (1952).)
 Membrane voltage is in absolute mV and has been reversed in polarity
  from the original HH convention and shifted to reflect a resting potential
  of -65 mV.
 Remember to set celsius=6.3 (or whatever) in your HOC file.
 See squid.hoc for an example of a simulation using this model.
 SW Jaslove  6 March, 1992
ENDCOMMENT
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
	(S) = (siemens)
}
 
? interface
NEURON {
        SUFFIX hh_Cp_scaled
        USEION na READ ena WRITE ina
        USEION k READ ek WRITE ik
        NONSPECIFIC_CURRENT il
        RANGE gnabar, gkbar, gl, el, gna, gk
        GLOBAL minf, ninf, mtau, ntau
}
 
PARAMETER {
        gnabar = .0 (S/cm2)	<0,1e9>
        gkbar = 0 (S/cm2)	<0,1e9>
        gl = 0 (S/cm2)	<0,1e9>
        el = -54.3 (mV)
}
 
STATE {
        m n
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
        minf ninf
	mtau (ms) ntau (ms)
}
 
LOCAL mexp, nexp        
 
? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gna = gnabar*m*m*m
	ina = gna*(v - ena)
        gk = gkbar*n*n*n*n
	ik = gk*(v - ek)      
        il = gl*(v - el)
}
 
 
INITIAL {
	rates(v)
	m = minf
	n = ninf
}

? states
DERIVATIVE states {  
        rates(v)
        m' =  (minf-m)/mtau
        n' = (ninf-n)/ntau
}
 
LOCAL q10


? rates
PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
        LOCAL  alpha, beta, sum
        TABLE minf, mtau, ninf, ntau DEPEND celsius FROM -100 TO 100 WITH 200

UNITSOFF
        q10 = 3^((celsius - 23)/10)
                :"m" sodium activation system
        alpha = -.182 * vtrap(-(v+55),6)
        beta =  -.124 * vtrap((v+55),6)
        sum = alpha + beta
	mtau = 0.25/(q10*sum)
        minf = alpha/sum
                :"n" potassium activation system
        alpha = .01*vtrap(-(v+55),10) 
        beta = .125*exp(-(v+65)/80)
	sum = alpha + beta
        ntau = 1/(q10*sum)
        ninf = alpha/sum
}
 
FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = -y*(1 - x/y/2)
        }else{
                vtrap = x/(1 - exp(x/y))
        }
}
 
UNITSON
