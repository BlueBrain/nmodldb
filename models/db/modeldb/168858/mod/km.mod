TITLE Potasium M type current for RD Traub, J Neurophysiol 89:909-921, 2003

COMMENT

	Implemented by Maciej Lazarewicz 2003 (mlazarew@seas.upenn.edu)

ENDCOMMENT

INDEPENDENT { t FROM 0 TO 1 WITH 1 (ms) }

UNITS { 
	(mV) = (millivolt) 
	(mA) = (milliamp) 
}
 
NEURON { 
	SUFFIX km
	USEION k READ ek WRITE ik
	RANGE gbar, ik, vshift
}

PARAMETER { 
	gbar = 0.0 	(mho/cm2)
	v ek 		(mV)  
        vshift = 0      (mV)
}
 
ASSIGNED { 
	ik 		(mA/cm2) 
	alpha beta	(/ms)
}
 
STATE {
	m
}

BREAKPOINT { 
	SOLVE states METHOD cnexp
	ik = gbar * m * ( v - ek ) 
}
 
INITIAL { 
	rates(v) 
	m = alpha / ( alpha + beta )
	m = 0
}
 
DERIVATIVE states { 
	rates(v) 
	m' = alpha * ( 1 - m ) - beta * m 
}

UNITSOFF 

PROCEDURE rates(v) { 
	alpha = 0.02 / ( 1 + exp( ( -v -vshift - 20 ) / 5 ) )
	beta  = 0.01 * exp( ( -v - vshift - 43 ) / 18 )
}

UNITSON
