TITLE CaL channel

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (S) = (siemens)
}

NEURON {
    THREADSAFE
    SUFFIX CaL
    USEION ca READ eca WRITE ica
    RANGE gmax, g, i
    GLOBAL vhalf, tmax
}

PARAMETER {
    gmax = 0.001 (S/cm2)
    vhalf = -15 (mV)
	vn2= -10 (mV)
	s1 = 3  (mV)
	s2 = 5  (mV)
	tmax = 0.5 (ms)
	tmin = 0.05 (ms)
} 

ASSIGNED {
    v (mV)
    eca (mV)
    
    minf
    tau (ms)
    ica (mA/cm2)
    i (mA/cm2)
	g (S/cm2)
}

STATE {
    m
}

INITIAL {
    settables(v)
    m = minf
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    g = gmax*m
    ica  = g*(v-eca)
    i=ica
}


DERIVATIVE states { 
    settables(v)      
    m' = (minf - m)/tau
}


PROCEDURE settables(v (mV)) {
    TABLE minf, tau DEPEND vhalf, tmax
          FROM -100 TO 50 WITH 1500

	: steady-state activation of ICaL in mV
    minf = 1/(1 + exp((vhalf-v)/s1))

    : steady-state ICaL time constant
	tau = tmin + tmax/(1 + exp((vn2-v)/s2))
}

