
NEURON {
POINT_PROCESS square
	RANGE del,dur,locx,locy,i,Vdependent
	RANGE g,gmax,m,e
	GLOBAL gama,n  
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) 	= (nanoamp)
	(mV)	= (millivolt)
	(nS) 	= (nanomho)
}
PARAMETER {
	del=50		(ms)
	dur=50		(ms)
	gmax=1		(nA)
	n=0.25 		(/mM)
	gama=0.08 	(/mV)
	e = 0 		(millivolts)
	locx=0		:location x
	locy=0		:location y
	Vdependent=1		:1 - voltage dependent 0- voltage independent
}

ASSIGNED {
	v 			(mV)
	i 			(nA)
	g           (nS)
	:m			(nS)
}

STATE {
	m
}

BREAKPOINT {
	IF(at_time(del)){
		state_discontinuity( m, gmax)
	}
	IF(at_time(del+dur)){
		state_discontinuity( m, 0)
	}
	IF (Vdependent==0){
		g=m
	}ELSE{
		g=m/(1+n*exp(-gama*v) )
	}
	i = (1e-3)*g * (v - e)
	
}

INITIAL {
	m=0
}
