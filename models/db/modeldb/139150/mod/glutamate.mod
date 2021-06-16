COMMENT
//****************************//
// Created by Alon Polsky 	//
//    apmega@yahoo.com		//
//		2010			//
//****************************//
ENDCOMMENT

TITLE NMDA synapse with depression

NEURON {
	POINT_PROCESS glutamate
	NONSPECIFIC_CURRENT inmda,iampa
	RANGE e ,gampamax,gnmdamax,local_v,inmda,iampa
	RANGE del,Tspike,Nspike
	RANGE gnmda,gampa
	GLOBAL n, gama,tau1,tau2,tau_ampa
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
	gnmdamax=1	(nS)
	gampamax=1	(nS)

	e= 0.0	(mV)
	tau1=90	(ms)	
	tau2=5	(ms)	
	tau_ampa=2	(ms)	
	n=0.25 	(/mM)	
	gama=0.08 	(/mV) 
	dt 		(ms)
	v		(mV)
	del=30	(ms)
	Tspike=10	(ms)
	Nspike=1

}

ASSIGNED { 
	inmda		(nA)  
	iampa		(nA)  
	gnmda		(nS)
	ica 		(mA/cm2)
	local_v	(mV):local voltage

}
STATE {
	A 		(nS)
	B 		(nS)
	gampa 	(nS)
}

INITIAL {
      gnmda=0 
      gampa=0 
	A=0
	B=0
}    

BREAKPOINT {  
    
	LOCAL count
	SOLVE state METHOD cnexp
	FROM count=0 TO Nspike-1 {
		IF(at_time(count*Tspike+del)){
			state_discontinuity( A, A+ gnmdamax)
			state_discontinuity( B, B+ gnmdamax)
			state_discontinuity( gampa, gampa+ gampamax)
		}
	}

	gnmda=(A-B)/(1+n*exp(-gama*v) )
	inmda =(1e-3)*gnmda*(v-e)
	iampa= (1e-3)*gampa*(v- e):*(1-1/(1+exp(-0.3*v)))
	local_v=v
}

DERIVATIVE state {
	A'=-A/tau1
	B'=-B/tau2
	gampa'=-gampa/tau_ampa
}





