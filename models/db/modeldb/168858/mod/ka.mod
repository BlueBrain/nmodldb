TITLE Potasium Type A current for RD Traub, J Neurophysiol 89:909-921, 2003

COMMENT

	Implemented by Maciej Lazarewicz 2003 (mlazarew@seas.upenn.edu)

ENDCOMMENT

INDEPENDENT { t FROM 0 TO 1 WITH 1 (ms) }

UNITS { 
	(mV) = (millivolt) 
	(mA) = (milliamp) 
} 
NEURON { 
	SUFFIX ka
	USEION k READ ek WRITE ik
	RANGE gbar, ik, vshift
}
PARAMETER { 
	gbar = 0.0 	(mho/cm2)
	v ek 		(mV)  
        vshift = 0.0    (mV)
} 
ASSIGNED { 
	ik 		(mA/cm2) 
	minf hinf 	(1)
	mtau htau 	(ms) 
} 
STATE {
	m h
}
BREAKPOINT { 
	SOLVE states METHOD cnexp
	ik = gbar * m * m * m * m * h * ( v - ek ) 
} 
INITIAL { 
	rates(v) 
	m  = minf
	m  = 0
	h  = hinf
} 
DERIVATIVE states { 
	rates(v) 
	m' = ( minf - m ) / mtau 
	h' = ( hinf - h ) / htau
}

UNITSOFF 

PROCEDURE rates(V (mV)) { 

	minf  = 1 / ( 1 + exp( ( - v -vshift- 60 ) / 8.5 ) )
	mtau = 0.185 + 0.5 / ( exp( ( v + vshift + 35.8 ) / 19.7 ) + exp( ( - v -vshift - 79.7 ) / 12.7 ) )
	hinf  = 1 / ( 1 + exp( ( v + vshift + 78 ) / 6 ) )
	if( v + vshift < -63 ) {
		htau = 0.5 / ( exp( ( v + vshift + 46 ) / 5 ) + exp( ( - v - vshift - 238 ) / 37.5 ) )
	}else{
		htau = 9.5
	}
}

UNITSON
