INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX cadyn
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
	depth	= 1	(um)		: depth of shell
	taur	= 100	(ms)		
	cainf	= 1e-4 (mM)
}

STATE {
	cai		(mM) 
}

INITIAL {
	cai = cainf
}

ASSIGNED {
	ica		(mA/cm2)
	drive_channel	(mM/ms)
}
	
BREAKPOINT {
	SOLVE state METHOD cnexp
}

DERIVATIVE state { 
	drive_channel =  - (10000)  * ica / (FARADAY * depth)

	cai' = drive_channel + (cainf-cai)/taur
}

