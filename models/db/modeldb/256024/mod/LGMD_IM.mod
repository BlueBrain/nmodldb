TITLE IM channel for LGMD


UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (S) = (siemens)
}

: gmax and g are range variables (i.e., can change in different compartments
NEURON {
    THREADSAFE
    : note - every variable accessible in NEURON will have the suffix _M
    SUFFIX M

    USEION k READ ek WRITE ik

    : these variables can be accessed as compartment.rangevar_M
    RANGE vhalf, gmax, tau, g

    : these will be accessed as taumax_M, vhalf_M, s1_M, s2_M  
    GLOBAL taumax, taumin, s1, s2
}

PARAMETER {
    gmax= 0.0003 (S/cm2)

    vhalf = -50 (mV)	: half activation
    s1 = 15		(mV)	: steepness of activation
	s2 = -20	(mV)
	taumax = 80 (ms)
	taumin = 15 (ms)
	aop = 0		(1)	< 0, 1 >	: voltage independent conductance
}

ASSIGNED { 
    v	(mV)
    ek	(mV)
    
    ik	(mA/cm2)
    ninf	(1)
    tau (ms)
    g 	(S/cm2)
}

STATE {
    n
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    g  = gmax*(n*(1-aop)+aop)
    ik  = g*(v-ek)
}

: calls the function settables below, then 
: set the steady state value of IM activation 
INITIAL {
    settables(vhalf-v)
    n = ninf
}

DERIVATIVE states { 
    settables(vhalf-v)      
    n' = (ninf - n)/tau
}


PROCEDURE settables(df (mV)) {
	: tabled relative to vhalf so one table can be used across sections with differing vhalf
    TABLE ninf, tau DEPEND s1, s2, taumax, taumin
          FROM -100 TO 100 WITH 600

    : steady-state activation of IM in mV
    ninf = 1/(1+(exp((df)/s1)))

    : steady-state IM time constant
    :tau = 4*taumax/(1+exp((df)/s2) )* ninf
    tau = 4*(taumax-taumin)/(1+exp((df)/s2))*ninf+taumin
}


