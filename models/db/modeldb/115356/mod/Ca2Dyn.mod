: adapted from My1stNRN
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX ca2dyn
	USEION ca READ ica
	USEION ca2 READ ca2i WRITE ca2i VALENCE 2.0
	RANGE depth,ca2inf,taur
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
	depth	= 1	(um)		: depth of shell
	taur	= 1000	(ms)		
	ca2inf	= 1e-4 (mM)
}

STATE {
	ca2i		(mM) 
}

INITIAL {
	ca2i = ca2inf
}

ASSIGNED {
	ica		(mA/cm2)
	drive_channel	(mM/ms)
}
	
BREAKPOINT {
	SOLVE state METHOD cnexp
}

DERIVATIVE state { 
	drive_channel =  - (10000)  * ica / (83.3333 * FARADAY * depth)

	ca2i' = drive_channel + (ca2inf-ca2i)/taur
}

