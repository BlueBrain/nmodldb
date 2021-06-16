TITLE K-Cl cotransporter KCC2

INDEPENDENT {t FROM 0 TO 1 WITH 10 (ms)}

NEURON {
	SUFFIX kcc2
	USEION k READ ko, ki WRITE ik
	USEION cl READ clo, cli WRITE icl VALENCE -1
	RANGE ik, icl
	GLOBAL U, S, Vi
}

UNITS {
	(mV)	= (millivolt)
	(molar) = (1/liter)
	(mM)	= (millimolar)
	(um)	= (micron)
	(mA)	= (milliamp)
	(mol)	= (1)
	FARADAY	= 96485.309 (coul/mole)
	PI	= (pi) (1)
}

PARAMETER {
	S = 5654.87 (um2)
  	Vi = 8913.48 (um3)
	U = 0.0003    (mM/ms)
}

ASSIGNED {
	ik		(mA/cm2)
	icl		(mA/cm2)
	ko		(mM)
	ki		(mM)
  	clo   (mM)
  	cli   (mM)
}

BREAKPOINT {
	LOCAL rate
	rate = pumprate(ki,ko,cli,clo)
	ik =  rate
	icl = -rate
}

FUNCTION pumprate(ki,ko,cli,clo) {
	pumprate = U*log((ki*cli)/(ko*clo))*(FARADAY*Vi/(S*1e4))
}
