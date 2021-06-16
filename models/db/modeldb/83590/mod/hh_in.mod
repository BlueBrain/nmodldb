TITLE hh_in.mod   squid sodium (slow inactivating), potassium, and leak channels
 
COMMENT

This is the original Hodgkin-Huxley treatment for the set of sodium,  potassium, and
leakage channels found in the squid giant axon membrane, incorporating the slow-cumulative
sodium inactivation, as reported by Miles et al., 2005.

PARAMETERS

 gnabar = .12   (S/cm2)    : maximal conductance for the sodium current.
 gkbar  = .036  (S/cm2)    : maximal conductance for the potassium current.
 gl     = .0003 (S/cm2)    : leak-current conductance.
 el     = -54.  (mV)       : reversal potential of the leak-current.
 a      = 1     (1)        : utility variable to switch on (i.e. a==1) and off (i.e. a==0) the slow-adaptation of sodium-current.

REFERENCES

Hodgkin, A.L., Huxley, A.F. (1952). A quantitative description of membrane current and its application conduction and excitation in nerve" J.Physiol. (Lond.) 117:500-544.
Miles, G.B., Dai, Y., and Brownstone, R.M. (2005). Mechanisms underlying the early phase of spike frequency adaptation in mouse spinal motoneurones. J Physiol 566.2 (2005) pp 519-532.
Arsiero, M., Luescher, H.-R., Lundstrom, B.N., and Giugliano, M. (2007). The Impact of Input Fluctuations on the Frequency-Current Relationships of Layer 5 Pyramidal Neurons in the Rat Medial Prefrontal Cortex. sumbitted.
  
AUTHORS

Michele Giugliano & Brian N. Lundstrom, Okinawa, June 5th 2006, and Lausanne Jan 5th 2007.
Modified from the original "hh.hoc" (SW Jaslove  6 March, 1992), provided with each standard distribution of NEURON.

ENDCOMMENT
 
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
         (S) = (siemens)
}
 
? interface
NEURON {
        SUFFIX hhin
        USEION na READ ena WRITE ina
        USEION k READ ek WRITE ik
        NONSPECIFIC_CURRENT il
        RANGE gnabar, gkbar, gl, el, gna, gk, a
        GLOBAL minf, hinf, ninf, sinf, mtau, htau, stau, ntau
}
 
PARAMETER {
        gnabar = .12   (S/cm2)    <0,1e9>
        gkbar  = .036  (S/cm2)    <0,1e9>
        gl     = .0003 (S/cm2)    <0,1e9>
        el     = -54.  (mV)
        a      = 1     (1)
}
 
STATE {
        m h n s
}
 
ASSIGNED {
        v       (mV)
        celsius (degC)
        ena     (mV)
        ek      (mV)

        gna (S/cm2)
        gk  (S/cm2)
        ina (mA/cm2)
        ik  (mA/cm2)
        il  (mA/cm2)
        minf hinf sinf ninf
    mtau (ms) htau (ms) stau (ms) ntau (ms)
}
 
:LOCAL mexp, hexp, sexp, nexp        
 
? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gna = gnabar*m*m*m*h*s
        ina = gna*(v - ena)
        gk = gkbar*n*n*n*n
        ik = gk*(v - ek)      
        il = gl*(v - el)
}
 
 
INITIAL {
    rates(v)
    m = minf
    h = hinf
    s = sinf
    n = ninf
}

? states
DERIVATIVE states {  
        rates(v)
        m' =  (minf-m)/mtau
        h' =  (hinf-h)/htau
        s' =  (sinf-s)/stau
        n' =  (ninf-n)/ntau
}
 
LOCAL q10


? rates
PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
        LOCAL  alpha, beta, sum
        TABLE minf, mtau, hinf, sinf, htau, ninf, stau, ntau DEPEND celsius FROM -100 TO 100 WITH 200


UNITSOFF
        q10 = 3^((celsius - 6.3)/10)
                :"m" sodium activation system
        alpha = .1 * vtrap(-(v+40),10) 
        beta =  4 * exp(-(v+65)/18)
        sum = alpha + beta
    mtau = 1/(q10*sum)
        minf = alpha/sum
                :"h" sodium inactivation system
        alpha = .07 * exp(-(v+65)/20)
        beta = 1 / (exp(-(v+35)/10) + 1)
        sum = alpha + beta
    htau = 1/(q10*sum)
        hinf = alpha/sum
                :"s" sodium inactivation system - according to Miles et al., 2005.
        alpha = 0.0077 / (1. + exp( (47.+v)/9. ) )
        beta  = 0.0077 / (1. + exp( -(47.+v)/9.  ) )
        sum = alpha + beta
    stau = 1/(sum)
        sinf = (1-a) + a*alpha/sum
                :"n" potassium activation system
        alpha = .01*vtrap(-(v+55),10) 
        beta = .125*exp(-(v+65)/80)
    sum = alpha + beta
        ntau = 1/(q10*sum)
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
