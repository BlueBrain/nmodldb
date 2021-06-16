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
    RANGE gmax, g, taus, tauh
    GLOBAL vhalf, tmax, thh, hmin
}

PARAMETER {
: all values can be adjusted in hoc files
    gmax = 0.001 (S/cm2)
    vhalf = -40 (mV)
    thh = -64	(mV)
    vn2= -50	(mV)
	s1 = 7.5	(mV)
	s2 = -8.0	(mV)
	tmax = 80	(ms)
	tmin = 5	(ms)
	kh = -5.4	(mV)
	hb = -100	(mV)
	hmin = 30	(ms)
	tadj = 20	(ms)
	hs = -14	(mV)
	nic = 0.0	(1)
	sp=1		(1)
}

ASSIGNED {
    v (mV)
    eca (mV)
    
    sinf (1)
    hinf (1)
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
    g = gmax*s^sp*h
    ica  = g*(v-eca)
}


DERIVATIVE states { 
    settables(v)
	s' = (sinf - s)/taus
	h' = (hinf - h)/tauh
}


PROCEDURE settables(v (mV)) {
    TABLE sinf, taus, hinf, tauh DEPEND vhalf, tmax, tmin, thh, hmin, tadj
          FROM -115 TO 50 WITH 660

	: steady-state activation of ICaT in mV
    sinf = 1/(1 + exp((vhalf-v)/s1))

    : steady-state CaT activation time constant
	taus = tmin + 4*(tmax-tmin)/(1 + exp((vn2-v)/s2))*sinf
	: add check to ensure positive tau
	
	: steady-state inactivation of ICaT in mV
	hinf = nic + (1-nic)/(1+exp((thh-v)/kh))

    : steady-state CaT inactivation time constant
	tauh = hmin + hinf*exp((hb-v)/hs)*tadj
	
}

