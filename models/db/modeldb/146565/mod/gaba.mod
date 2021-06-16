COMMENT
//****************************//
// Created by Alon Polsky 	//
//    apmega@yahoo.com		//
//		2010			//
//****************************//
based on Sun et al 2006
ENDCOMMENT

TITLE GABAA synapse 

NEURON {
	POINT_PROCESS gaba
	
	NONSPECIFIC_CURRENT i
	RANGE e ,gmax,local_v,i
	RANGE del,Tspike,Nspike
	RANGE taudgaba,dgaba,decaygaba
	RANGE R,D
	GLOBAL risetime,decaytime
}

PARAMETER {
	gmax=.5	(nS)
	e= -70.0	(mV)
	risetime=1	(ms)	:2
	decaytime=20(ms)	:40

	v		(mV)
	del=30	(ms)
	Tspike=10	(ms)
	Nspike=1
	taudgaba=200	(ms)
	decaygaba=0.6
	}

ASSIGNED { 
	i		(nA)  
	local_v	(mV):local voltage
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
}    

BREAKPOINT {  
    
	LOCAL count
	SOLVE state METHOD cnexp
	FROM count=0 TO Nspike-1 {
		IF(at_time(count*Tspike+del)){
			state_discontinuity( R, R+ gmax*dgaba)
			state_discontinuity( D, D+ gmax*dgaba)
			state_discontinuity( dgaba, dgaba* decaygaba)
		}
	}
	i= (D-R)*(1e-3)* (v- e)

	:i= (1e-3)*ggaba* (v- e)
	local_v=v
}

DERIVATIVE state {
	R'=-R/risetime
	D'=-D/decaytime
	dgaba'=(1-dgaba)/taudgaba

}





