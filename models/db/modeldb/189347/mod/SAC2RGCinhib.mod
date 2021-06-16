
:GABA release dependent on presynaptic voltage
NEURON {
POINT_PROCESS SACinhib
	RANGE Vpre ,Vinf
	GLOBAL tau ,Vtau ,e,maxves,gsingle ,newves
	RANGE release,numves,g,s_inf,t1,i,g
	RANGE locx,locy,local_v
	NONSPECIFIC_CURRENT i
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
	Vtau=30	(/ms)		:STIMULUS DEPOLARIZATION RATE	
							:postsynaptic
	gsingle=0.2	(nS)
	tau=10		(ms)
	e = -65 	(mV)
	locx=0		:location x
	locy=0		:location y

}

ASSIGNED {
	:presynaptic
	Vinf 		(mV)
	s_inf
	t1
	numves
	release
	:postsynaptic
	v 			(mV)
	i 			(nA)
	local_v		(mV)
}

STATE {
	g	 		(nS)
	Vpre 		(mV)
}
 
BREAKPOINT {
	SOLVE state METHOD euler
	if (t>t1){										:EVERY 1 MS
		releasefunc(Vpre)
		t1=t1+1
	}
	i = (1e-3)*g * (v - e)
	local_v=v
}

INITIAL {
	:presynaptic
	s_inf=0
	release=0
	numves=maxves
	t1=0
	Vinf=0
	Vpre=Vinf
	:postsynaptic
	g =0
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
		state_discontinuity( g, g+ release*gsingle)
	}
	addves=0							:REPLINISHMENT
	FROM rand=0 TO maxves-numves-1 {
		if (scop_random()<newves){addves=addves+1}
	}
	numves=numves+addves
}
DERIVATIVE state {
	g'=-g/tau
	Vpre'=(-Vpre+Vinf)/Vtau
}