TITLE IM channel for LGMD


UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (S) = (siemens)
}

: gmax and g are range variables (i.e., can change in different compartments
: while e is global
NEURON {
    THREADSAFE
    : note - every variable accessible in NEURON will be having the suffix _M
    SUFFIX M

    USEION k READ ek WRITE ik

    : these variables will be accessed as compartment.rangevar_M
    : note: to make the channel constant available add the following
    : to the next line: vhalf, s1, s2, tau_max
    RANGE gmax, g

    : this will be accessed as taumax_M, vhalf_M, s1_M, s2_M  
    GLOBAL taumax, taumin, vhalf, s1, s2
}

PARAMETER {
: all values can be adjusted in hoc files
    gmax= 0.001 (S/cm2)
    vhalf = -51	(mV)
    v2 = -60	(mV)
    s1 = 10.5	(mV)
    s2 = -11 	(mV)
    taumax = 60 (ms)
    taumin = 2	(ms)
}

ASSIGNED { 
    v (mV)
    ek (mV)
    
    ik (mA/cm2)
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
    ik = g*(v-ek)
}

: calls the function settables below, then 
: set the initial steady state value 
INITIAL {
    settables(v)
    n = ninf
}

DERIVATIVE states { 
    settables(v)      
    n' = (ninf - n)/tau
}


PROCEDURE settables(v (mV)) {
    :LOCAL vhalf, s1, s2, taumax
    :local variables take units of right hand side, see below

    TABLE ninf, tau DEPEND vhalf, s1, s2, taumax, taumin
          FROM -100 TO 50 WITH 750

    : steady-state activation of IM in mV
    ninf = 1/(1+(exp((vhalf-v)/s1)))

    : steady-state IM time constant
    : slope in mV and time constant in ms
	:tau = 4*(taumax-taumin)/(1+exp((v2-v)/s2))*ninf+taumin
	tau = 2*taumax/( exp((vhalf-v)/s2) + exp((vhalf-v)/s1) ) + taumin
}


