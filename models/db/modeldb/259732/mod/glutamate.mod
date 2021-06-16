COMMENT
//******************************************//
// Created by Alon Poleg-Polsky 			//
//    alon.poleg-polsky@ucdenver.edu		//
//		2018								//
//******************************************//
ENDCOMMENT

TITLE Glutamatergic synapse with network activation

NEURON {
	POINT_PROCESS glutamate
	NONSPECIFIC_CURRENT inmda,iampa
	RANGE e ,gAMPAmax,gNMDAmax,inmda,iampa

	RANGE gnmda,gampa,dend,pos,locx,locy,local_v
	RANGE stim,tt

	GLOBAL n, gama,tau_ampa,Pr
	GLOBAL tau1,tau2
	GLOBAL Voff,Vset

}

UNITS {
	(nA) 	= (nanoamp)
	(mV)	= (millivolt)
	(nS) 	= (nanomho)
	(mM)    = (milli/liter)
        F	= 96480 (coul)
        R       = 8.314 (volt-coul/degC)
 	PI = (pi) (1)
	(mA) = (milliamp)
	(um) = (micron)

}

PARAMETER {
	gNMDAmax=0.5	(nS)
	gAMPAmax=0.5	(nS)
	e= 0.0	(mV)
	tau1=50	(ms)	
	tau2=2	(ms)	
	tau_ampa=1	(ms)	
	n=0.25 	(/mM)	
	gama=0.08 	(/mV) 
	dt (ms)
	v		(mV)
	dend=0
	pos=0
	locx=0
	locy=0
	Pr=.8
	Voff      =0		:0 - voltage dependent 1- voltage independent
	Vset      =-60		:set voltage when voltage independent		
	stim=0	
}

ASSIGNED {
	inmda		(nA)  
	iampa		(nA)  
	gnmda		(nS)
	local_v
	tt
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
	stim=0
	tt=0
}    

BREAKPOINT {  
	if((stim==1)&&(tt<t)){
		if(scop_random()<=Pr){
			state_discontinuity( A, A+ gNMDAmax)
			state_discontinuity( B, B+ gNMDAmax)
			state_discontinuity( gampa, gampa+ gAMPAmax)
			tt=t+2
		}
	}
	SOLVE state METHOD cnexp
	local_v  =v*(1-Voff)+Vset*Voff	:temp voltage
	gnmda    =(A-B)/(1+n*exp(-gama*local_v) )
	inmda =(1e-3)*gnmda* (v-e)
	iampa= (1e-3)*gampa* (v- e)
}

DERIVATIVE state {
	A'=-A/tau1
	B'=-B/tau2
	gampa'=-gampa/tau_ampa
}
