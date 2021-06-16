COMMENT
//****************************//
// Created by Alon Polsky 	//
//    apmega@yahoo.com		//
//		2010			//
//****************************//
based on Sun et al 2006
ENDCOMMENT
TITLE GABAA synapse activated by the network
NEURON {
	POINT_PROCESS gaba_net
	NONSPECIFIC_CURRENT i
    RANGE xloc,yloc,tag1,tag2
	RANGE gmax,local_v,i
	RANGE taudgaba,dgaba,ggaba
	RANGE R,D
	GLOBAL risetime,decaytime ,e ,decaygaba
}
PARAMETER {
	gmax=.5	(nS)
	e= -60.0	(mV)
	risetime=1	(ms)	:2
	decaytime=20(ms)	:40

	v		(mV)
	taudgaba=200	(ms)
	decaygaba=0.8
	xloc=0
	yloc=0
	tag1=0
	tag2=0
}
ASSIGNED {
	i		(nA)  
 	local_v	(mV):local voltage
	ggaba
}

STATE {
	dgaba
	R
	D
}

INITIAL {
      dgaba=1 
	R=0
	D=0
	ggaba=0
}
BREAKPOINT {
	SOLVE state METHOD cnexp
	ggaba=D-R
	i=(1e-3)*(D-R)*(v-e)
	local_v=v
}
NET_RECEIVE(weight) {
	state_discontinuity( R, R+ gmax*dgaba)
	state_discontinuity( D, D+ gmax*dgaba)
	state_discontinuity( dgaba, dgaba* decaygaba)
}
DERIVATIVE state {
	R'=-R/risetime
	D'=-D/decaytime
	dgaba'=(1-dgaba)/taudgaba

}