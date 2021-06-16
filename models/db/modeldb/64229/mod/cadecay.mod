: $Id: cadecay.mod,v 1.1.1.1 2005/12/15 15:16:39 hines Exp $
TITLE First order calcium decay
: Paul Bush 3.31.92  No warranties expressed or implied.

INDEPENDENT {t FROM 0 TO 1 WITH 10 (ms)}

NEURON {
	SUFFIX cadecay
	USEION ca READ ica WRITE cai
	RANGE taucaremov
	GLOBAL cainit
}

UNITS {
	(molar) = (1/liter)
	(mM)	= (millimolar)
	(um)	= (micron)
	(mA)	= (milliamp)
:	FARADAY = 96520	(coul)
	FDY	= 9.6520	( 10000 coulomb)
}


PARAMETER {
	taucaremov	= 20	(ms)
	diam 		= 1	(um)
	ica			(mA/cm2)
	cainit		= 5e-5	(mM)
}

STATE {
	cai		(mM) 
}

BREAKPOINT {
	SOLVE state METHOD cnexp
}

DERIVATIVE state {
	cai' = ((cainit - cai)/taucaremov) + (-ica * 4/(diam*FDY))
: Note the 4 should actually be a 2 to take account of the double charge on 
: the calcium ion
}

INITIAL {
    cai = cainit
}
