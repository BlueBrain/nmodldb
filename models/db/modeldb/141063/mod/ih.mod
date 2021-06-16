TITLE Simplyfied model of h-current
: Written by Patricio Orio

NEURON { 
	SUFFIX ih
	NONSPECIFIC_CURRENT i
	RANGE  gbar
}

UNITS { 
	(mV) = (millivolt) 
	(mA) = (milliamp)
}

PARAMETER { 
	gbar = 0.0003 	(mho/cm2)
	tau = 100			(ms)
	z = -0.14		(/mV)
	V0 = -89		(mV)
	eh = -30   (mV)
} 

ASSIGNED { 
	g       	(mho/cm2)
	celsius		(degC)
	i 		(mA/cm2) 
	v 			(mV)
	ainf
	rho
	phi
}
 
STATE { a }

INITIAL { 
	a = 1 / (1+exp(-z*(v-V0)))
}

BREAKPOINT { 
	SOLVE kin METHOD cnexp
	rho = 1.3 ^((celsius - 25)/10 (degC) )
	i = rho*gbar*a*(v - eh)
}

DERIVATIVE kin {
	phi = 3^((celsius - 25)/10 (degC) )  
	ainf = 1 / (1+exp(-z*(v-V0)))
	a'=phi * (ainf - a)/tau
}
