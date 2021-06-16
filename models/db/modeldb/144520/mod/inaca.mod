TITLE sodium calcium exchange
COMMENT
	modified From DiFrancesco & Noble 1985 Phil Trans R Soc Lond 307:353-398 
    modified for Neuron by FE GANNIER
	francois.gannier@univ-tours.fr (University of TOURS)
ENDCOMMENT
INCLUDE "Unit.inc"
INCLUDE "Volume.inc"
NEURON {
	SUFFIX inaca
	USEION ca READ cai, cao WRITE ica
	USEION na READ nai, nao WRITE ina
	RANGE ina, ik, incx, Kncx
}

PARAMETER {
	Kncx = 0.02	(nA)
	nncx = 3
	dncx = 0.001
	gamma = 0.5
}

ASSIGNED {
	v (mV)
	celsius (degC) : 37
	ica (mA/cm2)
	ina (mA/cm2)
	incx (mA/cm2)
	cao (mM)
	nai (mM)
	cai (mM)
	nao (mM)
}

LOCAL RT
INITIAL {
	RT = (1000)*R*(273.15+celsius)
}

BREAKPOINT { 
	incx = (1e-06)* Kncx/S*(exp(gamma*(nncx-2)*v*F/RT)*(nai*1(/mM))^nncx*cao - exp((gamma-1)*(nncx-2)*v*F/RT)*(nao*1(/mM))^nncx*cai)/((1 + dncx*(cai*(nao*1(/mM))^nncx+cao*(nai*1(/mM))^nncx))*(1 + cai/0.0069(mM)))
	ina = 3*incx
	ica = -2*incx
}
