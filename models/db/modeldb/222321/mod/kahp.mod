TITLE kahp :Calcium-dependent K (afterhyperpolarization) current

NEURON { 
	SUFFIX kahp
	USEION k READ ek WRITE ik
	USEION ca READ cai
	RANGE gkahp, ik, m, qk
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
	gkahp = 0.0 	(mho/cm2) 
}
 
ASSIGNED { 
	ik      (mA/cm2) 
	alpha   (/ms)
	beta	(/ms)
	v	(mV)	
	ek 	(mV)
	diam	(um)
	cai	(mM)
}
 
STATE {	m qk }

BREAKPOINT { 
	SOLVE states METHOD cnexp
	ik = gkahp*m*(v-ek)
}

INITIAL { 
	rates(cai)
	m = alpha/(alpha+beta)
	ik = gkahp*m*(v-ek)
	qk = 0
}
 
DERIVATIVE states { 
	rates(cai)
	m' = alpha*(1-m)-beta*m
	qk' = (-ik*diam*PI*(1e4)/FARADAY)/(diam*diam*PI/4)
}

UNITSOFF

PROCEDURE rates(chi (mM)) { 

	if(cai<=1) {
		alpha = cai*30
	}else{
		alpha = 30
	}
	beta = 1.0
}

UNITSON
