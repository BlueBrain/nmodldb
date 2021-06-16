COMMENT
Internal calcium concentration calculated in response to a calcium current calculated elsewhere,
and a single-exponetial extruding pump calculated here. Valid for small distances from the membrane,
where the suface to volume ratio of a cylinder can be approximated as 1/depth .
This code was modified from Destexhe,Babloyantz & Sejnowski (Biophys. J. 65(1993): 1538-1552 --Equ. 7 ) 
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX cad
	USEION ca READ ica, cai WRITE cai
	RANGE depth,cainf,taur
}

UNITS {
	(molar) = (1/liter)			: moles do not appear in units
	(mM)	= (millimolar)
	(um)	= (micron)
	(mA)	= (milliamp)
	(msM)	= (ms mM)
}

CONSTANT {
	FARADAY = 96489		(coul)		: moles do not appear in units
}

PARAMETER {
	depth	= .1	(um)		: depth of shell
	taur	= 1.5	(ms)		: remove first-order decay
	cainf	= 0.0001 (mM)
}

STATE {
	cai		(mM) 
}

INITIAL {
	cai = 0.0001
}

ASSIGNED {
	ica		(mA/cm2)
	drive_channel	(mM/ms)
	drive_pump	(mM/ms)
}
	
BREAKPOINT {
	SOLVE state METHOD derivimplicit
}

DERIVATIVE state { 

: Surface to volume ratio of a cylinder is approximated as 1/depth for very small values of depth
: More precicely the ratio is   2 / ( (depth)^2 / (diam/2 - (2 * depth)) )

	drive_channel =  - (10000) * ica / (2 * FARADAY * depth)

	if (drive_channel <= 0.) { drive_channel = 0. }	: cannot pump below resting level
	cai' = drive_channel + (cainf-cai)/taur
}

