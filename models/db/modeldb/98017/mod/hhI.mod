TITLE hhI.mod   interneuron sodium, potassium, and leak channels
 
COMMENT

 This file is based on the original hh.mod file (see original comment
 below). It was modified to match the model that was used in the
 simulations of Wang and Buzsaki (1996, J. Neurosci. 16).

 Author: Fredrik Edin, 2003
 Address: freedin@nada.kth.se

 Original comment:
 ***************************************************************************
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
 ***************************************************************************

 changes:

  - m is substituted by its steady state value: m_inf - see 'BREAKPOINT'
  {as a result mtau is not needed, 'minf' is removed from
  GLOBAL declaration and 'm' is included in the RANGE var list
  otherwise it will be handled as a GLOBAL var and will not be
  evaluated separately for the 'sections'; for 'h' an 'n' this 
  is not a problem}

  - for h and n alpha and beta values are multiplied by 5 
  (see factor "Phi" in the W&B model)

  - temp: set to 6.3 Celsius, alpha and beta values are set/manipulated
  directly to simulate characteristic firing pattern

  ***************************************************************************
ENDCOMMENT
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
	(S) = (siemens)
}

? interface

NEURON {
        SUFFIX hhI
	USEION na READ ena WRITE ina
	USEION k READ ek WRITE ik
        NONSPECIFIC_CURRENT il
        RANGE gnabar, gna, gkbar, gk, gl, el, ek, ena, ina, ik, il
	GLOBAL hinf, ninf, htau, ntau, speed
}
 
PARAMETER {
	speed 	= 5 	(1)     	<0,1e9>
        gnabar	= .035 	(mho/cm2)	<0,1e9>
	gkbar 	= .009 	(mho/cm2)	<0,1e9>
        gl	= .0001 (mho/cm2)	<0,1e9>
        el 	= -65 	(mV)
}
  
STATE {
        h n
}
 
ASSIGNED {
        v (mV)
	celsius (degC)
	ena 	(mV)
	gna 	(mho/cm2)
        ina 	(mA/cm2) 
	ek	(mV)
	gk 	(mho/cm2) 
        ik 	(mA/cm2)
        il 	(mA/cm2)
        minf hinf ninf
	htau (ms) ntau (ms)
}
 
LOCAL mexp, hexp, nexp        
 
? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
	gna = gnabar*minf*minf*minf*h
	ina = gna*(v - ena)
        gk = gkbar*n*n*n*n
	ik = gk*(v - ek)      
        il = gl*(v - el)
}


INITIAL {
	rates(v)
	h = hinf
	n = ninf
}


? states
DERIVATIVE states {
        rates(v)
        h' = speed * (hinf-h)/htau
        n' = speed * (ninf-n)/ntau
}
 
LOCAL q10


? rates
PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
                          :Call once from HOC to initialize inf at resting v.
		      
        LOCAL  alpha, beta, sum
        TABLE minf, hinf, ninf, htau, ntau DEPEND celsius FROM -100 TO 100 WITH 200

UNITSOFF
        q10 = 3^((celsius - 6.3)/10)

        :"m" sodium activation system
        alpha = .1 * vtrap(-(v+35),10)
        beta =  4 * exp(-(v+60)/18)
        sum = alpha + beta
        minf = alpha/sum

        :"h" sodium inactivation system
        alpha = .07 * exp(-(v+58)/20)
        beta = 1 / (exp(-(v+28)/10) + 1)
        sum = alpha + beta
	htau = 1/(q10*sum)
        hinf = alpha/sum

        :"n" potassium activation system
        alpha = .01*vtrap(-(v+34),10) 
        beta = .125*exp(-(v+44)/80)
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
