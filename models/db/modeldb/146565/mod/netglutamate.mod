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
	RANGE Voff,Vset
	RANGE e ,gampamax,gnmdamax,local_v,inmda,iampa,local_ca
	RANGE decayampa,decaynmda,dampa,dnmda
	RANGE gnmda,gampa,xloc,yloc,tag1,tag2
	GLOBAL n, gama,tau_ampa,taudampa,taudnmda,tau1,tau2
	GLOBAL icaconst
	USEION ca READ cai WRITE ica
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
	gnmdamax=1	(nS)
	gampamax=1	(nS)
	icaconst =0.1:1e-6
	e= 0.0	(mV)
	tau1=50	(ms)	
	tau2=2	(ms)	
	tau_ampa=2	(ms)	
	n=0.25 	(/mM)	
	gama=0.08 	(/mV) 
	dt 		(ms)
	v		(mV)
:	del=30	(ms)
:	Tspike=10	(ms)
:	Nspike=1
	Voff=0		:0 - voltage dependent 1- voltage independent
	Vset=-60		:set voltage when voltage independent
	decayampa=.5
	decaynmda=.5
	taudampa=200	(ms):tau decay
	taudnmda=200	(ms):tau decay

	xloc=0
	yloc=0
	tag1=0
	tag2=0
}

ASSIGNED { 
	inmda		(nA)  
	iampa		(nA)  
	gnmda		(nS)
	local_v	(mV):local voltage
	local_ca	:local calcium,
	ica			(nA)
	cai

}
STATE {
	A 		(nS)
	B 		(nS)
	gampa 	(nS)
	dampa
	dnmda
}


INITIAL {
      gnmda=0 
      gampa=0 
	A=0
	B=0
	dampa=1
	dnmda=1
	ica=0
}    

BREAKPOINT {  
    
	LOCAL count
	SOLVE state METHOD cnexp
	local_v=v*(1-Voff)+Vset*Voff
	gnmda=(A-B)/(1+n*exp(-gama*local_v) )
	inmda =(1e-3)*gnmda*(v-e)
	iampa= (1e-3)*gampa*(v- e)
	local_v=v
	local_ca=cai
	ica=inmda*0.1/(PI*diam)*icaconst
	inmda=inmda*.9

}
NET_RECEIVE(weight) {
	state_discontinuity( A, A+ gnmdamax*(dnmda))
	state_discontinuity( B, B+ gnmdamax*(dnmda))
	state_discontinuity( gampa, gampa+ gampamax*dampa)
	state_discontinuity( dampa, dampa* decayampa)
	state_discontinuity( dnmda, dnmda* decaynmda)
}
DERIVATIVE state {
	A'=-A/tau1
	B'=-B/tau2
	gampa'=-gampa/tau_ampa
	dampa'=(1-dampa)/taudampa
	dnmda'=(1-dnmda)/taudnmda
}





