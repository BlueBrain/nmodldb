TITLE 	Potassium Current 
COMMENT
	Author: Ronald van Elburg 
	Taken from  model for fast spiking neuron in J. Tegner, A. Compte, 
	X.J. Wang, J. Neurosci. 22(20): 9053-9062, 2002
	
 	Modifications:
 	ID	Date		Authors				Email					Description
	M_001
ENDCOMMENT



UNITS {
: (Abbreviation)= (Unit)  
	(mV)        = (millivolt) 
	(mA)        = (milliamp) 
	(pS) 		= (picosiemens)
	(um) 		= (micron)
: Abbreviation 	= (Constant) (Unit)
} 

NEURON { 
	SUFFIX KPyr
	USEION k READ ek WRITE ik
	RANGE gbar, ik
	GLOBAL  v_table_min, v_table_max,phin
}

PARAMETER { 
:	Parameter	=Initial Value 	(Units)			Description
	gbar        =  	1           (pS/um2)
	v 	                        (mV)  
	ek 		                    (mV)  
	phin        =   1	        (1)
	
: Table settings
	v_table_min 		= -120			(mV)
	v_table_max 		= 100			(mV)
} 

ASSIGNED { 
:   Parameter   Units		  Description
	ik		        (mA/cm2) 
	nalpha          (1/ms)
	nbeta           (1/ms)
} 

STATE {
	n                 (1)
}

BREAKPOINT { 
	settables(v) 
	SOLVE states METHOD cnexp
	ik =(1e-4)*gbar * n * n * n * n * (v - ek)
} 

INITIAL { 
	settables(v) 
	n=nalpha/(nalpha+nbeta)
} 

DERIVATIVE states { 
	settables(v) 
	n' =phin* ( nalpha*(1-n) -nbeta*n)  
}


UNITSOFF
 
PROCEDURE settables(v (mV)) { 
	TABLE nalpha, nbeta FROM 	v_table_min  TO v_table_max WITH 961
	nalpha  = 0.032*vtrap(-(v+52),0.2)
	nbeta   = 0.5*exp(-0.025*(57+v))
	}

UNITSON


FUNCTION vtrap(x, k) {
  if (fabs(x) < 1e-6) {
    vtrap = 1/(k * exp(k*x))
  } else {
    vtrap = x / (exp(k*x) - 1)
  }
}
