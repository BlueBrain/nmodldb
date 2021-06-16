TITLE CaS channel

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (S) = (siemens)
}

NEURON {
    THREADSAFE
    SUFFIX CaS
    USEION ca READ eca WRITE ica
    RANGE gmax, g, taum, tauh
    GLOBAL vhalfm, tmax, vhalfh, hmin
}

PARAMETER {
: all values can be adjusted in hoc files
    gmax = 0.001 (S/cm2)
    vhalfm =-34	(mV)
    vhalfh =-33	(mV)
    vn2= -55	(mV)
	s1 = 10.0	(mV)
	s2 = -8.0	(mV)
	tmax = 25	(ms)
	tmin = 1	(ms)
	kh = -5.5	(mV)
	hb = -120	(mV)
	hmin = 5	(ms)
	tadj = 10	(ms)
	hs = -40	(mV)
	nic = 0.0	(1)
	mp=2		(1)
}

ASSIGNED {
    v (mV)
    eca (mV)
    
    minf (1)
    hinf (1)
    taum (ms)
    tauh (ms)
    ica (mA/cm2)
	g (S/cm2)
}

STATE {
    m
    h
}

INITIAL {
    settables(v)
    m = minf
    h = hinf
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    g = gmax*m^mp*h
    ica  = g*(v-eca)
}


DERIVATIVE states { 
    settables(v)
	m' = (minf - m)/taum
	h' = (hinf - h)/tauh
}


PROCEDURE settables(v (mV)) {
    TABLE minf, taum, hinf, tauh DEPEND vhalfm, tmax, tmin, vhalfh, hmin, tadj
          FROM -100 TO 50 WITH 600

	: steady-state activation of ICaT in mV
    minf = 1/(1 + exp((vhalfm-v)/s1))

    : steady-state CaT activation time constant
	taum = tmin + 4*(tmax-tmin)/(1 + exp((vn2-v)/s2))*minf
	: add check to ensure positive tau
	
	: steady-state inactivation of ICaT in mV
	hinf = nic + (1-nic)/(1+exp((vhalfh-v)/kh))

    : steady-state CaT inactivation time constant
	tauh = hmin + hinf*exp((hb-v)/hs)*tadj
	
}

