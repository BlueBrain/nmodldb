TITLE decay of submembrane calcium concentration
:
: Internal calcium concentration due to calcium currents and pump.
: Modified from Destexhe et al. 1994
: Differential equations.
:
: This file contains the following mechanisms:
: Simple first-order decay or buffering:
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
: Originally written by Alain Destexhe, Salk Institute, Nov 12, 1992
: Modified by Albert Gidon, Humboldt-Universitat zu Berlin 18/09/2019
:

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
	depth	= 1	(um)		: depth of shell
	taur	= 1e10	(ms)		: remove first-order decay
	cainf	= 2.4e-4 (mM)
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

	drive_channel =  - (10000) * ica / (2 * FARADAY * depth)

	if (drive_channel <= 0.) { drive_channel = 0. }	: cannot pump inward

	cai' = drive_channel +  (cainf-cai)/taur
}







