: Calcium handler using parameters from Table 1 in
: Wilson and Callaway, 2000, J Neurophysiol

NEURON {
	SUFFIX CaHandler2
	USEION ca READ cai, ica WRITE cai
	GLOBAL Pmax, beta, Dapp, mult
	RANGE cai_prime
	POINTER d2cai, cal
}

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
	(mM) = (milli/liter)
	(um) = (micrometer)
}

PARAMETER {
	Pmax = 2	(um/ms)
	beta = 0.001
	mult = 1
	z = 2
	F = 9.64846e4	(coulomb)
	Dapp = .02	(um2/ms)
}

ASSIGNED {
	diam		(um)
	cai			(mM)
	ica			(mA/cm2)
	icalcium	(mA/m2)
	dt			(ms)
	cai_prime	(mM/ms)
	d2cai		(mM/um2)
	cal			(mM)
}


BREAKPOINT {
	rate()
	caistep()
}


INITIAL {
	rate()
	cai = cal
}

PROCEDURE rate() {
	icalcium = ica * 1e4 (cm2/m2)
	cai_prime = (Dapp*d2cai) -((icalcium*4*beta*mult)/(z*F*diam)) - ((Pmax*4*beta*cal)/diam)
}

PROCEDURE caistep() {

	cai = cal + cai_prime*dt
}