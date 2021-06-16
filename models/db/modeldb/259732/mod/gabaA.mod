COMMENT
//**********************************//
// Created by Alon Poleg-Polsky 	//
// alon.poleg-polsky@ucdenver.edu	//
// 2018								//
//**********************************//
ENDCOMMENT

TITLE GABAergic synapse with network activation

NEURON {
	POINT_PROCESS gabaA
	NONSPECIFIC_CURRENT i
	RANGE e ,gmax
	RANGE g,dend,pos,locx,locy,local_v
	GLOBAL tau
	RANGE stim,tt,Pr
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
	gmax=0.5	(nS)
	e= -70.0	(mV)	
	tau=7	(ms)	
	dt (ms)
	v		(mV)
	dend=0
	pos=0
	locx=0
	locy=0
	stim=0	
	Pr=.8		
}

ASSIGNED {
	i		(nA)  
	local_v
	tt
}
STATE {
	g 	(nS)
}

INITIAL {
    g=0 
	i=0
	stim=0
	tt=0	
}    

BREAKPOINT {  
	if((stim==1)&&(tt<t)){
		if(scop_random()<=Pr){
			state_discontinuity( g, g+ gmax)
			tt=t+2
		}
	}
	SOLVE state METHOD cnexp
	i= (1e-3)*g* (v- e)
	local_v=v
}

DERIVATIVE state {
	g'=-g/tau
}
