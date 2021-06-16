TITLE  H-current that uses Na ions

:  Updated to use CVode - Carl Gold  08/10/03
:  Updated by Maria Markaki  03/12/03

NEURON {
	SUFFIX h
	USEION na READ ena WRITE ina      
        RANGE  gbar,vhalf, K, taun, ninf,ina
:	NONSPECIFIC_CURRENT i
}

UNITS {
	(um) = (micrometer)
	(mA) = (milliamp)
	(uA) = (microamp)
	(mV) = (millivolt)
}

PARAMETER {              : parameters that can be entered when function is called in cell-setup
        eh     = -10   (mV)
	K      = 8.5   (mV)
	gbar   = 0     (mho/cm2)  : initialize conductance to zero
	vhalf  = -90   (mV)       : half potential
}	

ASSIGNED {             : parameters needed to solve DE
	v              (mV)
        ena            (mV)
	ina            (mA/cm2)
	ninf
	taun           (ms)
}

STATE {                : the unknown parameters to be solved in the DEs
	n
} 


INITIAL {               : initialize the following parameter using states()
	states()	
	n = ninf
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	ina = gbar*n*(v-eh)            
}

DERIVATIVE states {
	rates(v)
        n' = (ninf - n)/taun
}

PROCEDURE rates(v (mV)) {  
 
 	if (v > -30) {
	   taun = 1(ms)
	} else {
           taun = 2(ms)*(1/(exp((v+145)/-17.5(mV))+exp((v+16.8)/16.5(mV))) + 5) :h activation tau

	}  
         ninf = 1 - (1 / (1 + exp((vhalf - v)/K)))                  :steady state value
}



