COMMENT

	Calcium activated K channel from Av-Ron and Vidal, 1999
	Implemented by C. Weaver, 2003

ENDCOMMENT

UNITS {
	(molar) = (1/liter)
}

UNITS {
	(mV) =	(millivolt)
	(mA) =	(milliamp)
	(mM) =	(millimolar)
}


NEURON {
	SUFFIX kca
	USEION ca READ cai
	USEION k READ ek WRITE ik
	RANGE gbar,gkca
	RANGE tot
}

PARAMETER {
	celsius		(degC)
	v		(mV)
	gbar=.001	(mho/cm2)	: Maximum Permeability
	: cai = 5.e-5	(mM)

	Kd=0.5		(mM)
}

ASSIGNED {
	ik		(mA/cm2)
	tot (mA/cm2)
      gkca          (mho/cm2)
	ek		(mV)
	cai		(mM)
}

INITIAL {
	gkca = gbar*cai/(Kd+cai)
: printf("kca ik=%g\n",ik)
}

BREAKPOINT {
	gkca = gbar*cai/(Kd+cai)
	tot = gkca*(v - ek)
	ik = gkca*(v - ek)
}


