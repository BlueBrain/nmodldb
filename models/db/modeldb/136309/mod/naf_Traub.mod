TITLE Sodium transient current for RD Traub et al 2003, 2005

COMMENT

	Implemented by Maciej Lazarewicz 2003 (mlazarew@seas.upenn.edu)
	fastNashift init to 0 and removed from arg modification Tom Morse 3/8/2006
	(for Traub et al 2005)
ENDCOMMENT

INDEPENDENT { t FROM 0 TO 1 WITH 1 (ms) }

UNITS { 
	(mV) = (millivolt) 
	(mA) = (milliamp) 
} 
NEURON { 
	SUFFIX nafTraub
	USEION na READ ena WRITE ina
	RANGE gbar, ina,m, h, df, minf, mtau
}
PARAMETER { 
	: putting as written fastNa_shift = 0: orig -3.5 (mV)
	gbar = 0.0 	   (mho/cm2)
	v (mV) ena 		   (mV)  
} 
ASSIGNED { 
	ina 		   (mA/cm2) 
	minf hinf 	   (1)
	mtau (ms) htau 	   (ms)
	df	(mV)
} 
STATE {
	m h
}
BREAKPOINT { 
	SOLVE states METHOD cnexp
	ina = gbar * m * m * m * h * ( v - ena ) 
	df = v - ena
} 
INITIAL { 
	settables( v )
	m = minf
	m = 0
	h  = hinf
} 
DERIVATIVE states { 
	settables( v ) 
	m' = ( minf - m ) / mtau 
	h' = ( hinf - h ) / htau
}

UNITSOFF 

PROCEDURE settables(v1(mV)) {

	TABLE minf, hinf, mtau, htau  FROM -120 TO 40 WITH 641

	minf  = 1 / ( 1 + exp( (- v1 - 34.5 ) / 10 ) )
	if( v1 < -26.5 ) {
		mtau = 0.025 + 0.14 * exp( ( v1 + 26.5 ) / 10 )
	} else{
		mtau = 0.02 + 0.145 * exp( ( - v1 - 26.5 ) / 10 ) 
	}

	: hinf, and htau are shifted 3.5 mV comparing to the paper

	hinf  = 1 / ( 1 + exp( ( v1 + 59.4 ) / 10.7 ) )
	htau = 0.15 + 1.15 / ( 1 + exp( ( v1 + 33.5 ) / 15 ) )
}

UNITSON
