TITLE Potasium AHP type current for RD Traub, J Neurophysiol 89:909-921, 2003

COMMENT

     Written by Aniruddha Yadav (anirddha.yadav@mssm.edu)

ENDCOMMENT

UNITS { 
	(mV) = (millivolt) 
	(mA) = (milliamp) 
} 

NEURON { 
	SUFFIX skahp
	USEION k READ ek WRITE ik
	USEION ca READ cai
	RANGE gbar, ik
}

PARAMETER { 
	gbar = 0.0 	(mho/cm2)
	v		(mV) 
	ek 		(mV)  
	cai		(1)
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
	rates( cai )
	m = alpha / ( alpha + beta )
	m = 0
}
 
DERIVATIVE states { 
	rates( cai )
	m' = alpha * ( 1 - m ) - beta * m 
}

UNITSOFF 

PROCEDURE rates(chii) { 

	if( cai < 100 ) {
		alpha = cai / 100000
	}else{
		alpha = 0.001
	}
	beta = 0.001
}

UNITSON
