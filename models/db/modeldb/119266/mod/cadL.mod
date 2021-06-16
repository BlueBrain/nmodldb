
NEURON {
	SUFFIX cadL
	USEION cal READ ical, cali WRITE cali  VALENCE 2	
	GLOBAL depth,calinf,taurl
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
	taurl	= 2.5	(ms)		: rate of calcium removal
	calinf	= 1e-4(mM)
}

STATE {
	cali		(mM) 
}


ASSIGNED {
	ical		(mA/cm2)
	drive_channel	(mM/ms)
}
	
BREAKPOINT {
	SOLVE state METHOD derivimplicit
}

DERIVATIVE state { 

	drive_channel =  - (10000) * ical / (2 * FARADAY * depth)
	if (drive_channel <= 0.) { drive_channel = 0.  }   : cannot pump inward 
         
	cali' = drive_channel + (calinf-cali)/taurl
}
