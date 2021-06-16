TITLE Ih
: Ih current from "Barnes and Hille (1989)"
:
: Rodrigo Publio 
: Revised - 20/01/2006

NEURON {
	SUFFIX ih
	NONSPECIFIC_CURRENT ih
	RANGE ghbar, eh
	GLOBAL ninf,taun
}


PARAMETER {
	gbar =  0.0311  	(mho/cm2)
	eh   = 	-17.7		(mV)            
	v 			(mV)
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	
} 

ASSIGNED {
	ih 		(mA/cm2)
	ninf	
	taun 		(ms)
}
 

STATE { n }

BREAKPOINT {
        SOLVE states METHOD cnexp
	ih = gbar*(1-(1+3*n)*(1-n)^3)*(v - eh)
} 

INITIAL {
	trates(v)
	n=ninf  
}

DERIVATIVE states {   
        trates(v)      
        n' = (ninf-n)/taun
}

PROCEDURE trates(vm) {  
        LOCAL  a, b

	a = alpha(vm,110,3,15)
	b = beta(vm,115,1.5,15)
	ninf = a/(a+b)
	taun = (1/(a+b))*100
}

FUNCTION alpha(v,vhalf,k,s) {
	
	        alpha = k / ((exp((v + vhalf)/s))+1)
	
}	

FUNCTION beta(v,vhalf,k,s) {
	
	        beta = k / ((exp(-(v + vhalf)/s))+1)
	
}





   