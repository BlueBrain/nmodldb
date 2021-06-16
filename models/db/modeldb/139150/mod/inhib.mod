COMMENT
//****************************//
// Created by Alon Polsky 	//
//    apmega@yahoo.com		//
//		2009			//
//****************************//
ENDCOMMENT

TITLE Shunting inhibition

NEURON {
	POINT_PROCESS shunt
	
	NONSPECIFIC_CURRENT ishunt
	RANGE e ,gmax,local_v,ishunt
	RANGE del,dur
	RANGE gshunt
}

UNITS {
	(nA) 	= (nanoamp)
	(mV)	= (millivolt)
	(nS) 	= (nanomho)
	(mM)    = (milli/liter)


}

PARAMETER {
	gmax=1	(nS)
	e= -70.0	(mV)
	dt (ms)
	v		(mV)
	del=0	(ms)
	dur=100 (ms)
}

ASSIGNED { 
	ishunt		(nA)  
	local_v	(mV):local voltage
}
STATE {
	gshunt

}

INITIAL {
      gshunt=0 

}    

BREAKPOINT {  
    
	IF(at_time(del)){
		state_discontinuity( gshunt, gmax)
	}
	IF(at_time(del+dur)){
		state_discontinuity( gshunt, 0)
	}	

	ishunt= (1e-3)*gshunt* (v- e)
	local_v=v
}







