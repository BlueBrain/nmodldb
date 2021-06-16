
TITLE Sodium transient current for dendritic channels showing little inactivation modified RD Traub, J Neurophysiol 89:909-921, 2003

COMMENT

	Implemented by Maciej Lazarewicz 2003 (mlazarew@seas.upenn.edu)

ENDCOMMENT

INDEPENDENT { t FROM 0 TO 1 WITH 1 (ms) }

UNITS { 
	(mV) = (millivolt) 
	(mA) = (milliamp) 
} 
NEURON { 
	SUFFIX na
	USEION na READ ena WRITE ina
	RANGE gbar, ina
}
PARAMETER { 
	fastNashift2 = -17 (mV) : changed from -3.5
	gbar = 0.0 	   (mho/cm2)
	v ena 		   (mV)  
} 
ASSIGNED { 
	ina 		   (mA/cm2) 
	minf2 hinf2 	   (1)
	mtau2 htau2 	   (ms) 
} 
STATE {
	m2 h2
}
BREAKPOINT { 
	SOLVE states METHOD cnexp
	ina = gbar * m2 * m2 * m2 * h2 * ( v - ena ) 
} 
INITIAL { 
	settables( v - fastNashift2 ) 
	m2 = minf2
	m2 = 0
	h2 = hinf2
} 
DERIVATIVE states { 
	settables( v ) 
	m2' = ( minf2 - m2 ) / mtau2 
	h2' = ( hinf2 - h2 ) / htau2
}

UNITSOFF 

PROCEDURE settables(v1(mV)) {

	TABLE minf2, hinf2, mtau2, htau2  FROM -120 TO 40 WITH 641

	minf2  = -0.02 + 1 / ( 1 + exp( ( - ( v1 + fastNashift2 ) - 55 ) / 8 ) )
	if( ( v1 + fastNashift2 ) < -40.0 ) {
		mtau2 = 0.025 + 0.14 * exp( ( ( v1 + fastNashift2 ) + 55 ) / 10 ) : original Tau = 10
 	} else{
		mtau2 = 0.01 + 0.45 * exp( ( - ( v1 + fastNashift2 ) - 55 ) / 10 ) : makes the cell need more current, original Tau = 10
	}

	: hinf, and htau are shifted 3.5 mV comparing to the paper

	hinf2  = 1 / ( 1 + exp( ( ( v1 + fastNashift2 * 0 ) + 62.9 ) / 10.7 ) ) : Original value 0.08
	htau2 = 0.15 + 1.15 / ( 1 + exp( ( ( v1 + fastNashift2 * 0)+ 37 ) / 10 ) ) : Original value 0.15
}

UNITSON



