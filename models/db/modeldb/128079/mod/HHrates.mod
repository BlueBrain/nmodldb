: Eight state kinetic scheme for HH sodium channel
: Modified from k3st.mod, chapter 9.9 (example 9.7)
: of the NEURON book
: 12 August 2008, Christoph Schmidt-Hieber

VERBATIM
#include <math.h>
ENDVERBATIM

NEURON {
    THREADSAFE
    SUFFIX HHrates
    USEION na READ ena WRITE ina
    GLOBAL vShift, maxrate, vShift_inact
    RANGE g, gbar, am0, am1, am2, bm0, bm1, ah0, ah1, bh0, bh1, bh2, hScale
}

UNITS { (mV) = (millivolt) }

: initialize parameters according to Engel & Jonas 2005
PARAMETER {
    gbar = 33     (millimho/cm2)
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
    
    vShift = 0        (mV) :shift to the right to account for Donnan potentials 
    vShift_inact = 0       (mV)  : global additional shift to the right for inactivation
    hScale = 1.0           :scales alpha_h and beta_h
    maxrate = 8.00e+03 (/ms) : limiting value for reaction rates
}

ASSIGNED {
    v    (mV)
    ena  (mV)
    g    (millimho/cm2)
    ina  (milliamp/cm2)
    am   (/ms)
    bm   (/ms)
    ah   (/ms)
    bh   (/ms)
}

STATE { c1 c2 c3 i1 i2 i3 i4 o }

BREAKPOINT {
    SOLVE kin METHOD sparse
    g = gbar*o
    ina = g*(v - ena)*(1e-3)
}

INITIAL { SOLVE kin STEADYSTATE sparse }

KINETIC kin {
    rates(v)
    ~ c1 <-> c2 (3*am, 1*bm)
    ~ c2 <-> c3 (2*am, 2*bm)
    ~ c3 <-> o  (1*am, 3*bm)
    ~ i1 <-> i2 (3*am, 1*bm)
    ~ i2 <-> i3 (2*am, 2*bm)
    ~ i3 <-> i4 (1*am, 3*bm)
    ~ i1 <-> c1 (ah, bh)
    ~ i2 <-> c2 (ah, bh)
    ~ i3 <-> c3 (ah, bh)
    ~ i4 <-> o  (ah, bh)
    CONSERVE c1 + c2 + c3 + i1 + i2 + i3 + i4 + o = 1
}

FUNCTION_TABLE tau1(v(mV)) (ms)
FUNCTION_TABLE tau2(v(mV)) (ms)

FUNCTION my_exp(x) {   :Non-overflowing exp
    if (x>700) {
	my_exp = pow( exp(1.0), x )
    } else {
	my_exp = exp(x)
    }
}

FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
    if (fabs(x/y) < 1e-6) {
        vtrap = y*(1 - x/y/2)
    }else{
        vtrap = x/( my_exp(x/y) - 1 )
    }
}

PROCEDURE rates(v(millivolt)) {
    LOCAL vS
    
    vS = v-vShift

    am = am0 * vtrap(-(vS-am1), am2)
    am = am*maxrate / (am+maxrate)
    bm = bm0 * my_exp(-vS/bm1)
    bm = bm*maxrate / (bm+maxrate)

    ah = hScale * ah0 * my_exp((-vS-vShift_inact)/ah1)
    ah = ah*maxrate / (ah+maxrate)
    bh = hScale * bh0 / (my_exp(-(vS-vShift_inact+bh1)/bh2) + 1)
    bh = bh*maxrate / (bh+maxrate)
}
