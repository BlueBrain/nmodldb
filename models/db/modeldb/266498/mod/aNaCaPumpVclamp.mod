TITLE sodium potassium pump
:  Adapted from Leo Medina's implementation from Lindblad et al Am J Physiol 1996 275:H1666

: Original model has been modified to assume constant nai

NEURON {
	SUFFIX aNaCaPump
	USEION ca READ cao, cai WRITE ica	
	USEION na READ nao, nai WRITE ina
	RANGE  inca, DFout, DFin, S, KNaCa, DNaCa
	GLOBAL dummy2, na_int : prevent vectorization for use with CVODE
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
  	(mM) = (millimolar)
	F = (faraday) (kilocoulombs)
	R = (k-mole) (joule/degC)
}

PARAMETER {
	KNaCa = 1.27324E-06 (mA/cm2/mM4)    <0,1e6>
	na_int = 8.9 (mM)
	Q10NaCa = 2.20
	r=3
	gamma=0.5
	DNaCa=0.0036 (/mM4)
}

ASSIGNED {
	celsius (degC)
	v (mV)
	cai (mM)
	cao (mM)
	ica (mA/cm2)
	ina (mA/cm2)
	nao (mM)
	nai (mM)
	inca (mA/cm2)
	dummy2
	S
	DFin (mM4)
	DFout (mM4)
	temp (degC)
}

BREAKPOINT {

	temp = celsius +273.15

	S=1.0+DNaCa*(cai*nao*nao*nao+cao*nai*nai*nai)
	
	DFin=nai*nai*nai*cao*exp(((r-2)*gamma*v*F)/(R*temp))
	
	DFout=nao*nao*nao*cai*exp(((r-2)*(gamma-1)*v*F)/(R*temp))

	inca=KNaCa*((DFin-DFout)/S)
	
		if (celsius >= 37) {
		inca=Q10NaCa*inca
	}
	
	ina = 3*inca
	ica = -2*inca
}