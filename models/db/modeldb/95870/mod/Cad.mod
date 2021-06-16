
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX Cad
	USEION Ca READ iCa, Cai WRITE Cai,Cao VALENCE 2	
        RANGE Ca, depth, Cainf, taur, entryF
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
	taur	= 200	(ms)		: rate of calcium removal
	Cainf	= 2e-3  (mM)		: 2uM
	Cai		(mM)
	Cao     = 2     (mM)
	entryF  = 1
}

STATE {
	Ca		(mM) 
}

INITIAL {
	Ca = Cainf
	Cao=2
	
}

ASSIGNED {
	iCa		(mA/cm2)
	drive_channel	(mM/ms)
}
	
BREAKPOINT {
	SOLVE state METHOD derivimplicit
}

DERIVATIVE state { 

	drive_channel =  - (10000) * iCa / (2 * FARADAY * depth)
	if (drive_channel <= 0.) { drive_channel = 0.  }   : cannot pump inward 
         
	:ca' = drive_channel + (cainf-ca)/taur
        :Ca' = drive_channel/18 + (Cainf -Ca)/taur*7
	Ca' = entryF*drive_channel/2 + (Cainf-Ca)/taur
	
        Cai = Ca
	Cao=2 :mM
}
