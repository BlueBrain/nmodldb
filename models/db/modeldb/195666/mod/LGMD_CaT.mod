TITLE CaT channel

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (S) = (siemens)
}

NEURON {
    THREADSAFE
    SUFFIX CaT
    USEION ca READ eca WRITE ica
    RANGE gmax, g
    GLOBAL vhalf, tmax, thh, hmin
}

PARAMETER {
: all values can be adjusted in hoc files
    gmax = 0.001 (S/cm2)
    vhalf = -56 (mV)
    thh = -71	(mV)
    vn2= -52	(mV)
	s1 = 8.0	(mV)
	s2 = -10	(mV)
	tmax = 4.4	(ms)
	tmin = 2.2	(ms)
	kh = -4.1	(mV)
	hb = -130	(mV)
	hmin = 9	(ms)
	tadj = 1	(ms)
	hs = -11.5	(mV)
	nic = 0.0	(1)
}

ASSIGNED {
    v (mV)
    eca (mV)
    
    sinf
    hinf
    taus (ms)
    tauh (ms)
    ica (mA/cm2)
	g (S/cm2)
}

STATE {
    s
    h
}

INITIAL {
    settables(v)
    s = sinf
    h = hinf
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    g = gmax*s^3*h
    ica  = g*(v-eca)
}


DERIVATIVE states { 
    settables(v)
	s' = (sinf - s)/taus
	h' = (hinf - h)/tauh
}


PROCEDURE settables(v (mV)) {
    TABLE sinf, taus, hinf, tauh DEPEND vhalf, tmax, thh, hmin
          FROM -150 TO 50 WITH 800

	: steady-state activation of ICaT in mV
    sinf = 1/(1 + exp((vhalf-v)/s1))

    : steady-state CaT activation time constant
	taus = tmin + 4*(tmax-tmin)/(1 + exp((vn2-v)/s2))
	: add check to ensure positive tau
	
	: steady-state inactivation of ICaT in mV
	hinf = nic + (1-nic)/(1+exp((thh-v)/kh))

    : steady-state CaT inactivation time constant
	tauh = hmin + hinf*exp((hb-v)/hs)*tadj
	
}

