TITLE hhmfb.mod   modified squid sodium, potassium, and leak channels

COMMENT
  This is a Hodgkin-Huxley treatment for the set of sodium, 
  potassium, and leakage channels found in the squid giant axon membrane, 
  adapted to channels in mossy fiber terminals
  Original: SW Jaslove  6 March, 1992
  Modified version: P Jonas, 10 March, 2004
  Modified version: C Schmidt-Hieber, 17 Dec 2007
  Modified version: C Schmidt-Hieber, 08 Sep 2008
  Modified version: C Schmidt-Hieber, 17 Apr 2010
  Revision history: 
  September 13, 2004: Revised fit - final version
  Dec 17, 2007: hScale to account for slower inactivation in the soma, CSH
  Sep 08, 2008: Make rate parameters range variables
  Apr 17, 2010: global inactivation shift; make vShift (Donnan) global
ENDCOMMENT

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
}
 
? interface
NEURON {
    THREADSAFE
    SUFFIX hhmfb
    USEION na READ ena WRITE ina
    USEION k READ ek WRITE ik
    NONSPECIFIC_CURRENT il
    RANGE gnabar, gkbar, gl, el, gna, gk, hScale, am0, am1, am2, bm0, bm1, ah0, ah1, bh0, bh1, bh2
    GLOBAL minf, hinf, ninf, mtau, htau, ntau, vShift, vShift_inact
}
 
PARAMETER {
    : alpha = 93.8285 * vtrap(-(v-105.023-vShift+vLeft), 17.7094)
    am0  = 9.38285e+1 (/ms)
    am1  = 1.05023e+2 (mV) : Note that this is used as a positive value here.
    am2  = 1.77094e+1 (mV)
    : beta =  0.168396 * exp(-(v-vShift+vLeft)/23.2707)
    bm0  = 1.68396e-1 (/ms)
    bm1  = 2.32707e+1 (mV)
    : alpha = hScale * .000353747 * exp(-(v-vShift)/18.706)
    ah0  = 3.53747e-4 (/ms)
    ah1  = 1.87060e+1 (mV)
    : beta = hScale * 6.62694 / (exp(-(v+17.6769-vShift)/13.3097) + 1)
    bh0  = 6.62694e+0 (/ms)
    bh1  = 1.76769e+1 (mV)
    bh2  = 1.33097e+1 (mV)
    
    gnabar = .12 (mho/cm2)	<0,1e9>
    gkbar = .036 (mho/cm2)	<0,1e9>
    gl = .0003 (mho/cm2)	      <0,1e9>
    el = -80.0 (mV)
    vShift = 12 (mV) :shift to the right to account for Donnan potentials 
    vShift_inact = 0 (mV) :global inactivation shift to align with 8-state model
    hScale = 1   : account for slower inactivation in the soma
}
 
STATE {
        m h n
}
 
ASSIGNED {
    v (mV)
    ena (mV)
    ek (mV)
    gna (mho/cm2)
    gk (mho/cm2)
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
    rates(v)
    m = minf
    h = hinf
    n = ninf
}

? states
DERIVATIVE states {  
    rates(v)
    m' =  (minf-m)/mtau
    h' = (hinf-h)/htau
    n' = (ninf-n)/ntau
}
 

? rates
PROCEDURE rates(v(mV)) {  
    :Computes rate and other constants at current v.
    :Call once from HOC to initialize inf at resting v.
    LOCAL  alpha, beta, sum, vS
    TABLE minf, mtau, hinf, htau, ninf, ntau FROM -100 TO 100 WITH 200

UNITSOFF
    vS = v-vShift
    : "m" sodium activation system
    : am1 is positive here (by contrast to the original E&J model)
    alpha = am0 * vtrap(-(vS-am1), am2)
    beta =  bm0 * exp(-vS/bm1)
    sum = alpha + beta
    mtau = 1/sum
    minf = alpha/sum

    :"h" sodium inactivation system
    alpha = hScale * ah0 * exp(-(vS-vShift_inact)/ah1)
    beta = hScale * bh0 / (exp(-((vS-vShift_inact)+bh1)/bh2) + 1)
    sum = alpha + beta
    htau = 1/sum
    hinf = alpha/sum
    
    :"n" potassium activation system
    alpha = .01*vtrap(-(v+55),10) 
    beta = .125*exp(-(v+65)/80)
    sum = alpha + beta
    ntau = 1/sum
    ninf = alpha/sum
}
 
FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
    if (fabs(x/y) < 1e-6) {
        vtrap = y*(1 - x/y/2)
    } else {
        vtrap = x/(exp(x/y) - 1)
    }
}
 
UNITSON
