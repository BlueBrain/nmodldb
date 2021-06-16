TITLE Na-K Pump

INDEPENDENT {t FROM 0 TO 1 WITH 10 (ms)}

NEURON {
	SUFFIX nakpump
	USEION k READ ko WRITE ik
	USEION na READ nai WRITE ina
	RANGE ik, ina, km_k, km_na, imax, qna, qk
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
	R 	= (k-mole)(joule/degC)
}

PARAMETER {
	celsius		(degC)
	km_k = 2	(mM) 
	km_na = 10	(mM)
	imax = 0.013    (mA/cm2)
}

STATE { qna qk }

ASSIGNED {
	ik		(mA/cm2)
	ina		(mA/cm2)
	ko		(mM)
	nai		(mM)
	diam		(um2)
	L		(um)
}

BREAKPOINT {
	SOLVE integrate METHOD sparse
}

INITIAL {
	qna=0
	qk=0
	ik = -2*imax*flux(nai,ko)
	ina = ik*-3/2
}

KINETIC integrate {
	ik = -2*imax*flux(nai,ko)
	ina = ik*-3/2

	COMPARTMENT diam*diam*PI/4 { qna qk }
	~ qna << (-ina*PI*diam*(1e4)/FARADAY)
	~ qk  <<  ( -ik*PI*diam*(1e4)/FARADAY)
}

FUNCTION flux(na,k) {
	flux = (1/((1+km_k/k)*(1+km_k/k)))*(1/((1+km_na/na)*(1+km_na/na)*(1+km_na/na)))
}
