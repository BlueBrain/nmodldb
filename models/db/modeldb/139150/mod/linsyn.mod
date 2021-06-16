COMMENT
//****************************//
// Created by Alon Polsky 	//
//    apmega@yahoo.com		//
//		2007			//
//****************************//
ENDCOMMENT

TITLE linear synapse; can be (AMPA/GABAa) synapse 

NEURON {
	POINT_PROCESS linsyn
	
	NONSPECIFIC_CURRENT i
	RANGE e ,gmax,local_v
	RANGE del,Tspike,Nspike
	RANGE gsyn
	RANGE tau_syn
}

UNITS {
	(nA) 	= (nanoamp)
	(mV)	= (millivolt)
	(nS) 	= (nanomho)
	(mM)    = (milli/liter)
        F	= 96480 (coul)
        R       = 8.314 (volt-coul/degC)

}

PARAMETER {
	gmax=1	(nS)
	e= 0.0	(mV)
	tau_syn=2	(ms)	

	dt (ms)
	v		(mV)
	del=30	(ms)
	Tspike=10	(ms)
	Nspike=1

}

ASSIGNED { 
	i		(nA)  
	local_v	(mV):local voltage
}
STATE {
	gsyn

}

INITIAL {
      gsyn=0 

}    

BREAKPOINT {  
    
	LOCAL count
	SOLVE state METHOD cnexp
	FROM count=0 TO Nspike-1 {
		IF(at_time(count*Tspike+del)){
			state_discontinuity( gsyn, gsyn+ gmax)
		}
	}

	i= (1e-3)*gsyn* (v- e)
	local_v=v
}

DERIVATIVE state {
	gsyn'=-gsyn/tau_syn
}





