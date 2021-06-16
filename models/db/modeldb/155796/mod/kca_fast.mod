TITLE C-type potasium current 
: Taken from RD Traub, J Neurophysiol 89:909-921, 2003
: Implemented by Maciej Lazarewicz 2003 (mlazarew@seas.upenn.edu)
: Adapted calcium dependence by Jordan Chambers 2012 (jordandchambers@gmail.com)

INDEPENDENT { t FROM 0 TO 1 WITH 1 (ms) }

UNITS { 
	(mV) = (millivolt) 
	(mA) = (milliamp) 
	(mM) = (milli/liter)
}
 
NEURON { 
	SUFFIX kca_fast
	USEION k READ ek WRITE ik
	USEION ca READ cai
	RANGE  gbar, ik, jikcaf, jk1, jk2
}

PARAMETER { 
	gbar = 3e-2 	(mho/cm2)
	v (mV) 
	ek 		(mV)  
	cai		(mM)
	cac = 6e-4 (mM)
	cas = 7.5e-5 (mM)
	v1 = 50
	v2 = 53.5
	s1 = 11
	s2 = 27
	eK=  -95 (mV)
	jikcaf		(mA/cm2)
	jk1
	jk2
	aspeed = 1
	bspeed = 1
} 

ASSIGNED { 
	ik 		(mA/cm2) 
	alpha (/ms) beta	(/ms)
}
 
STATE {
	m
}

BREAKPOINT { 
     SOLVE states METHOD cnexp
     jk1 = 1/(1+exp(((cac-cai)/cas)))
     jk2 = gbar*m*(v-ek)
     ik = gbar*m*(1/(1+exp(((cac-cai)/cas))))*(v-ek)
     jikcaf = ik
}
 
INITIAL { 
	settables(v) 
	m = alpha / ( alpha + beta )
	m = 0
}
 
DERIVATIVE states { 
	settables(v) 
	m' = alpha * ( 1 - m ) - beta * m 
}

UNITSOFF 

PROCEDURE settables(v(mV)) { 
	TABLE alpha, beta FROM -120 TO 40 WITH 641

	if( v <= -10.0 ) {
		alpha = 2 / 37.95 * ( exp( ((v + v1)/s1) - ((v + v2)/s2)))

		: Note that there is typo in the paper - missing minus sign in the front of 'v'
		beta  = 2 * exp((- v - v2)/s2) - alpha
	}else{
		alpha = 2 * exp(( - v - v2)/s2)
		beta  = 0
	}
	: speed-up of C kinetics here.
	alpha = alpha * aspeed
	beta  = beta  * bspeed
}

UNITSON
