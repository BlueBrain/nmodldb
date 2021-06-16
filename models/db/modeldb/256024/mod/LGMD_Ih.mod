TITLE Ih channel for LGMD
: this channel is identical to that implemented in Matlab code
: based on Richard's current and voltage clamp data (Sept 13)
: this is a second version were the state variables are ninf and tau 
: instead of nalpha and nbeta to allow for inspection of these variables 
: during simulations
: Modified on 04/08/14 to correct the apparent mix-up between tau and tau_max
: tau_max was RANGE and tau GLOBAL when it should be the opposite. 

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (S) = (siemens)
}

: gmax and g are range variables (i.e., can change in different compartments
: while e is global
NEURON {
    THREADSAFE
    : note - every variable accessible in NEURON will be having the suffix _h
    SUFFIX h

    NONSPECIFIC_CURRENT i

    : these variables will be accessed as compartment.rangevar_h
    : note: to make the channel constant available add the following
    : to the next line: vhalf, s1, s2, tau_max
    RANGE gmax, tau, g, taumax, i

    : this will be accessed as e_h, taumax_h, vhalf_h, s1_h, s2_h 
    GLOBAL e, taumin, vhalf, s1, s2
}

PARAMETER {
    gmax= 0.001 (S/cm2)
    e = -35 (mV)

    : note the following is only in the case we define these parameters as accessible to 
    : NEURON. Otherwise it is sufficient to initialize them in the procedure settables
    vhalf = -78 (mV)
    s1 = -13 (mV)
    s2 = 14 (mV)
    taumax = 1350 (1)
    taumin = 10 (ms)
}

ASSIGNED { 
    v (mV)

    i (mA/cm2)
    ninf
    tau (ms)
    g (S/cm2)
}

STATE {
    n
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    g  = gmax*n
    i  = gmax*n*(v-e)
}

: calls the function settables below, then 
: set the steady state value of Ih activation 
INITIAL {
    settables(v)
    n = ninf
}

DERIVATIVE states { 
    settables(v)      
    n' = (ninf - n)/(taumax*tau+taumin)
}


PROCEDURE settables(v (mV)) {
    :LOCAL vhalf, s1, s2, taumax
    :local variables take units of right hand side, see below

    TABLE ninf, tau DEPEND vhalf, s1, s2
          FROM -200 TO 50 WITH 750

    : steady-state activation of Ih in mV
    :vhalf = -77.8 (mV)
    :s1 = 13.8 (mV)
    ninf = 1/(1+(exp((vhalf-v)/s1)))

    : steady-state Ih time constant
    : slope in mV and time constant in ms
    :s2 = 19.7 (mV)
    :taumax = 1071.1 (ms)
    :tau = 2*taumax/( exp((v-vhalf)/s2) + exp((vhalf-v)/s2) )
    tau = (4 (ms))/(1+exp((vhalf-v)/s2))*ninf

}


