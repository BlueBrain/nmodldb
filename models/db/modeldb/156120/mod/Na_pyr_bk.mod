TITLE 	Sodium Current 
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
	SUFFIX NaPyr
	USEION na READ ena WRITE ina
	RANGE  gbar, ina, minf,malpha, mbeta
	GLOBAL  v_table_min, v_table_max
}

PARAMETER {
:	Parameter	=Initial Value 	(Units)			Description
	gbar =  	350     (pS/um2)
	v                   (mV)
	ena 		        (mV)  
	phih =  5	        (1)
	phim =	5           (1)
	
:   Table Settings
	v_table_min 		= -120			(mV)
	v_table_max 		= 100			(mV)
} 

ASSIGNED { 
:   Parameter   Units		  Description
	ina 		        (mA/cm2) 
	:minf 		        (1)			:Steady state activation approximation
	halpha              (1/ms)
	hbeta               (1/ms)
    malpha          	(1/ms)      :Dynamic activation
   	mbeta           	(1/ms)      :Dynamic (de)activation
} 


STATE {
	h                 (1)
	m                 (1)
}

BREAKPOINT { 
	settables(v) 
	SOLVE states METHOD cnexp
	:ina =(1e-4)*gbar * minf * minf * minf * h * (v - ena)  :Steady state activation approximation
    ina =(1e-4)*gbar * m * m* m * h * (v - ena)            :Dynamic activation
} 

INITIAL {
	settables(v) 
	h=halpha/(halpha+hbeta) 
	m=malpha/(malpha+mbeta)
} 

DERIVATIVE states { 
	h' =phih* ( halpha*(1-h) -hbeta*h) 
    m' =phim* ( malpha*(1-m) - mbeta*m)             :Dynamic activation
}


UNITSOFF


PROCEDURE settables(v (mV)) { 
	TABLE halpha, hbeta, malpha, mbeta FROM 	v_table_min  TO v_table_max WITH 961   :For dynamic activation add malpha, mbeta; For steady state activation approximation add minf
	
	halpha  = 0.128*(exp(-(v+50)/18))
	hbeta   = 4/(1+exp(-0.2*(v+27)))
	
	
	malpha = 0.32*vtrap(-(v+54),0.25)   
	mbeta =  0.28*vtrap((v+27),0.2)                      
	
	:minf  =  malpha/(malpha+mbeta) 					:Steady state activation approximation
}

UNITSON


FUNCTION vtrap(x, k) {
  if (fabs(x) < 1e-6) {
    vtrap = 1/(k * exp(k*x))
  } else {
    vtrap = x / (exp(k*x) - 1)
  }
}
