
NEURON {
	SUFFIX cadN
	USEION can READ ican, cani WRITE cani  VALENCE 2	
	GLOBAL depth,caninf,taurn
}

UNITS {
	(molar) = (1/liter)			: moles do not appear in units
	(mM)	= (millimolar)
	(um)	= (micron)
	(mA)	= (milliamp)
	(msM)	= (ms mM)
	FARADAY = (faraday) (coulomb)
}


PARAMETER {
	depth	= .1	(um)		: depth of shell
	taurn	= 1	(ms)		: rate of calcium removal
 	caninf	= 0.8e-4(mM)
}

STATE {
	cani		(mM) 
}


ASSIGNED {
	ican		(mA/cm2)
	drive_channel	(mM/ms)
}
	
BREAKPOINT {
	SOLVE state METHOD derivimplicit
}

DERIVATIVE state { 

	drive_channel =  - (10000) * ican / (2 * FARADAY * depth)
	if (drive_channel <= 0.) { drive_channel = 0.  }   : cannot pump inward 
         
	cani' = drive_channel + (caninf-cani)/taurn
}
