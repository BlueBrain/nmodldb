
:5403813100

:glutamate (AMPA+NMDA) release dependent on presynaptic voltage modeled here
NEURON {
POINT_PROCESS bipNMDA
	RANGE Vpre,Vdel,Vdur,Vamp,Vbase,locx,locy,local_v,i,g,release,numves :Vpre1,Vpre2,
	RANGE gAMPA,gNMDA,s_inf,t1,A,B,Vinf	:1,Vinf2
	GLOBAL maxves,newves,gAMPAsingle,gNMDAsingle,Vtau,Voff,Vset,VampK	:,Vnoiserate,Vnoiseamp,Vnoise,Vtau2
	GLOBAL icaconst ,gama,n,e,tauAMPA,tau1NMDA ,tau2NMDA
	NONSPECIFIC_CURRENT iAMPA,iNMDA
	:USEION ca WRITE ica
}

UNITS {
	(nA) 	= (nanoamp)
	(mV)	= (millivolt)
	(nS) 	= (nanomho)
}
PARAMETER {
							:presynaptic
	maxves=10				:TOTAL NUMBER OF VESICLES
	newves=0.01				:REPLENISHMENT RATE - VESICLES
	:Vnoise=0				:VOLTAGE VARIABILITY-INPUT NOISE
	Vdel=50 	(ms)		:START OF ACTIVATION
	Vdur=100	(ms)		:DURATION OF ACTIVATION
	Vamp=10		(mV)		:PEAK AMPLITUDE -STIMULUS
	Vbase=0		(mV)		:BACKGROUND ACTIVATION
	VampK=2					:SIZE OF THE INITIAL RESPONSE
	:Vnoiseamp=.1			:INTRINSIC NOISE
	:Vnoiserate=50			:INTRINSIC NOISE
	Vtau=30	(/ms)		:STIMULUS DEPOLARIZATION RATE
	:Vtau2=100	(/ms)		:STIMULUS HYPERPOLARIZATION RATE
							:postsynaptic
	gAMPAsingle=0.2	(nS)	:AMPA CONDUCTANCE
	gNMDAsingle=0.2	(nS)	:NMDA CONDUCTANCE
	tau1NMDA=50	(ms)		:DEACTIVATION
	tau2NMDA=2	(ms)		:ACTIVATION
	tauAMPA=2	(ms)		:DEACTIVATION
	n=0.25 		(/mM)		:NMDA VOLTAGE DEPENDENCE 
	gama=0.08 	(/mV)		:NMDA VOLTAGE DEPENDENCE 
	e = 0 		(mV)		:REVERSAL POTENTIAL
	locx=0					:location x
	locy=0					:location y
	icaconst =0.1			:CALCIUM FRACTION
	Voff=0					:0 - voltage dependent 1- voltage independent
	Vset=-60				:set voltage when voltage independent
}

ASSIGNED {
	:presynaptic
	Vinf 		(mV)
	:Vinf2 		(mV)	
	s_inf
	t1	
	numves
	release

	:postsynaptic
	v 			(mV)
	i 			(nA)
	g           (nS)
	iNMDA		(nA)
	iAMPA		(nA)
	gNMDA		(nS)
	local_v		(mV)
	:ica			(nA)
}

STATE {
	A
	B
	gAMPA 		(nS)
	:Vpre1 		(mV)
	:Vpre2 		(mV)
	Vpre		(mV)
}

BREAKPOINT {
	SOLVE state METHOD euler
	if (t>t1){										:EVERY 1 MS
		:Vpre=Vpre1-Vpre2+Vbase						:+normrand(0,Vnoise)	:PRESYNAPTIC VOLTAGE
		:if (Vpre<0){
		:	Vpre=0
		:}
:		if (scop_random()<1/Vnoiserate){			:ADDED INTRINSIC NOISE
:			Vinf1=Vinf1*(1+Vnoiseamp*(scop_random()-.5))
:		}
		releasefunc(Vpre)
		t1=t1+1
	}
	:IF(at_time(Vdel)){								:START STIMULUS 
	:	Vinf1=Vamp*VampK									:*(1+Vnoiseamp*(scop_random()-.5))
	:	Vinf2=Vamp*(VampK-1)
	:}
	:IF(at_time(Vdel+Vdur)){							:END STIMULUS
	:	Vinf1=Vbase									:*(1+Vnoiseamp*(scop_random()-.5))
	:	Vinf2=Vbase
	:}

	local_v=v*(1-Voff)+Vset*Voff					:VOLTAGE DEPENDENCE
	gNMDA=(A-B)/(1+n*exp(-gama*local_v) )
	iAMPA = (1e-3)*gAMPA * (v - e)
	iNMDA = (1e-3)*gNMDA * (v - e)
	i= iAMPA+iNMDA									:INDICATOR OF TOTAL CURRENT
	g=gNMDA+gAMPA									:INDICATOR OF TOTAL CONDUCTANCE

	:ica=iNMDA*icaconst	
	:iNMDA=iNMDA-ica
}

INITIAL {
	:presynaptic
	s_inf=0
	release=0
	numves=maxves
	t1=0
	Vinf=0	:Vbase
	Vpre=Vinf
	:Vinf2=Vbase
	:Vpre1=Vbase
	:Vpre2=Vbase

	:postsynaptic
	gAMPA=0
	gNMDA=0
	A=0
	B=0
}
 
FUNCTION releasefunc(vpre){
	LOCAL rand,addves
	s_inf=vpre/100
	release=0	
	FROM rand=0 TO numves-1 {			:GOES OVER ALL RRP
		if (scop_random()<s_inf){
			release=release+1
		}
	}
	if (release>0){						:RELEASE
		numves=numves-release
		if (numves<0){numves=0}
		state_discontinuity( gAMPA, gAMPA+ release*gAMPAsingle)
		state_discontinuity( A, A+ release*gNMDAsingle)
		state_discontinuity( B, B+ release*gNMDAsingle)
	}
	addves=0							:REPLINISHMENT
	FROM rand=0 TO maxves-numves-1 {
		if (scop_random()<newves){addves=addves+1}
	}
	numves=numves+addves
}
DERIVATIVE state {
	A'=-A/tau1NMDA
	B'=-B/tau2NMDA
	gAMPA'=-gAMPA/tauAMPA
	Vpre'=(-Vpre+Vinf)/Vtau
	:Vpre2'=(-Vpre2+Vinf2)/Vtau2

}