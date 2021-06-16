TITLE km :Slowly activated voltage-dependent K current

NEURON {
	SUFFIX km
	USEION k READ ek WRITE ik
	RANGE gkm, ik, alpha, beta
}

UNITS { 
	(molar) = (1/liter)
	(mV) =	(millivolt)
	(mA) =	(milliamp)
	(mM) =	(millimolar)
        FARADAY	= 96485.309 (coul/mole)
	PI	= (pi) (1) 
}

INDEPENDENT {t FROM 0 TO 1 WITH 100 (ms)}

PARAMETER { 
	gkm = 0.0 	(mho/cm2)
} 

ASSIGNED { 
	ik 	(mA/cm2) 
	alpha   (/ms)
        beta	(/ms)
	ek 	(mV)
	diam	(um)
	v       (mV)
	cai	(mM)
}
 
STATE {	m qk }

BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = gkm*m*(v-ek)
}

INITIAL {
	settables(v) 
	m = alpha/(alpha+beta)
	ik = gkm*m*(v-ek)
	qk = 0
}

DERIVATIVE states { 
	settables(v) 
	m' = alpha*(1-m)-beta*m
	qk' = (-ik*diam*PI*(1e4)/FARADAY)/(diam*diam*PI/4)
}

UNITSOFF 

PROCEDURE settables(v (mV)) { 
	TABLE alpha, beta FROM -120 TO 40 WITH 641

	alpha = m_inf(v)/tau_act(v)
	beta = (1-m_inf(v))/tau_act(v)
}

UNITSON

FUNCTION m_inf(v(mV)) {
	m_inf = 1/(1+exp(-(v+35)/5))
}

FUNCTION tau_act(v(mV)) {
	tau_act = 1000/(3.3*exp((v+35)/40)+exp(-(v+35)/20))
}
