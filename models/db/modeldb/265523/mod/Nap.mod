TITLE nap

NEURON {
	THREADSAFE
	SUFFIX nap
	USEION na READ ena WRITE ina
	RANGE  gbar, thegna, sh
	GLOBAL minf, mtau 
	:, hinf, mtau, htau
}

PARAMETER {
	gbar = .0006   	(mho/cm2)
	sh = 0  (mV)
	eNa = 55 	(mV)		
	ena		(mV)            
	celsius (degC)
	v 		(mV)
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

ASSIGNED {
	ina 		(mA/cm2)
	thegna		(mho/cm2)
	minf 	
	mtau (ms)
}
 

STATE { m }

UNITSOFF

BREAKPOINT {
    SOLVE states METHOD cnexp
	        		
	thegna =gbar*m       
	ina = thegna * (v - eNa)
	} 

INITIAL {
	mtau = 5
	minf = (1/(1+exp(-(v+51-sh)/5)))      	
	m=minf  
	
}

DERIVATIVE states {   
    
	mtau = 5
	minf = (1/(1+exp(-(v+51-sh)/5))) 	         	
	m' = (minf-m)/mtau
}



UNITSON
