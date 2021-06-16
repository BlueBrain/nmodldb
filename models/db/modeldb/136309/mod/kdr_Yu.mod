TITLE Potasium dr type current for Yu et al. 2008

COMMENT

	Implemented by Erin Munro

ENDCOMMENT

INDEPENDENT { t FROM 0 TO 1 WITH 1 (ms) }

UNITS { 
	(mV) = (millivolt) 
	(mA) = (milliamp) 
} 

NEURON { 
	SUFFIX kdrYu
	USEION k READ ek WRITE ik
	RANGE gbar, ik, m, am, bm, dfkdr
}

PARAMETER { 
	gbar = 0.0 	(mho/cm2)
	v 		(mV)  
	ek		(mV)
}
 
ASSIGNED { 
	ik 		(mA/cm2) 
	am 		(1/ms)
	bm 		(1/ms) 
	df		(mV)
}
 
STATE {
	m
}

BREAKPOINT { 
	SOLVE states METHOD cnexp
	df = v - ek
	ik = gbar * m * ( df ) 
}
 
INITIAL { 
	settables(v) 
	m = 0
}
 
DERIVATIVE states { 
	settables(v) 
	m' = am*(1-m) - bm*m
}

UNITSOFF 

PROCEDURE settables(v (mV)) { 
	TABLE am, bm FROM -120 TO 40 WITH 641

  am = 0.02*(v-25)/(1-exp(-(v-25)/9))
  bm = -0.002*(v-25)/(1-exp((v-25)/9))
}

UNITSON
