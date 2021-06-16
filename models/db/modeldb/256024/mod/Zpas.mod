TITLE passive membrane with capacitance and inductance

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
	(S) = (siemens)
	(uF) = (microfarad)
}

NEURON {
    THREADSAFE

	SUFFIX Zpas
	NONSPECIFIC_CURRENT i
	RANGE g0, e, C, L, pl
}

PARAMETER {
	g0 = 1e-4	(S/cm2)	<0,1e9>
	e = -65		(mV)
	C = 1.0		(uF/cm2)	<0,1e3>
	L = 2.5e3	(henry-cm2)	<0,1e9>
	pl = 0.2	(1)	<0,1>
}

ASSIGNED {
	v (mV)
	i (mA/cm2)
	tauC (ms)
	tauL (ms)
}

STATE {
	cv	(mV)
	lv	(mV)
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	
	tauL = L*g0*pl*1e3(ms/s)
    if (tauL<1e-2)	{ tauL = 1e-2 (ms) }
	tauC = C/g0*1e-3(ms/us)
    if (tauC<1e-2)	{ tauC = 1e-2 (ms) }
    
	i = g0*(1-pl)*(v-e) + g0*pl*(lv-e) + g0*(cv-v)
}

INITIAL {
	cv = v
	lv = v
}

DERIVATIVE states {
	lv' = (v - lv)/tauL
	cv' = (v - cv)/tauC
}
