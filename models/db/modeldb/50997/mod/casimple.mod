TITLE decay of submembrane calcium concentration
:
: 2. Simple first-order decay or buffering:
:
:       Cai + B <-> ...
:
:   which can be written as:
:
:       dCai/dt = (cainf - Cai) / taur
:
:   where cainf is the equilibrium intracellular calcium value (usually
:   in the range of 200-300 nM) and taur is the time constant of calcium 
:   removal.  The dynamics of submembranal calcium is usually thought to
:   be relatively fast, in the 1-10 millisecond range (see Blaustein, 
:   TINS, 11: 438, 1988).
:
: All variables are range variables
:
: Written by Alain Destexhe, Salk Institute, Nov 12, 1992
: Modified by TJ Velte, Univ of Minnesota, May 6, 1995

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX cas
	USEION ca READ ica, cai WRITE cai
	RANGE kd,cainf,taur
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
	diam (micron)
	taur	= 1.5	(ms)		: remove first-order decay
	cainf	= 0.0001 (mM)
	kd	= 0.0011	(mM)
}

STATE {
	cai		(mM) 
}

INITIAL {
	cai = kd
}

ASSIGNED {
	ica		(mA/cm2)
	drive_channel	(mM/ms)
}
	
BREAKPOINT {
	SOLVE state METHOD derivimplicit
					: printf("diam=%g\n", diam)
}

DERIVATIVE state { 

	drive_channel = - ( (2 * ica) / (FARADAY * diam) ) 

	if (drive_channel <= 0.) { drive_channel = 0. }	: cannot pump inward

	cai' = drive_channel - ( (cai - cainf) /  taur )
}

