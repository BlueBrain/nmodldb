COMMENT
//****************************//
// Created by Alon Polsky 	//
//    apmega@yahoo.com		//
//		2007			//
//****************************//
ENDCOMMENT

TITLE NMDA synapse with depression

NEURON {
	POINT_PROCESS ampa
	
	NONSPECIFIC_CURRENT iampa
	RANGE e ,gmax,ntar,local_v,iampa,gh
	RANGE del,Tspike,Nspike
	RANGE gampa
	GLOBAL tau_ampa
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
	tau_ampa=2	(ms)	

	dt (ms)
	v		(mV)
	del=30	(ms)
	Tspike=10	(ms)
	Nspike=1

}

ASSIGNED { 
	iampa		(nA)  
	local_v	(mV):local voltage
}
STATE {
	gampa

}

INITIAL {
      gampa=0 

}    

BREAKPOINT {  
    
	LOCAL count
	SOLVE state METHOD cnexp
	FROM count=0 TO Nspike-1 {
		IF(at_time(count*Tspike+del)){
			state_discontinuity( gampa, gampa+ gmax)
		}
	}

	iampa= (1e-3)*gampa* (v- e)
	local_v=v
}

DERIVATIVE state {
	gampa'=-gampa/tau_ampa
}





