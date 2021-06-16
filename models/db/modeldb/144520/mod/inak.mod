TITLE Na+/K+ pump current
COMMENT
	modified From DiFrancesco & Noble 1985 Phil Trans R Soc Lond 307:353-398 
    modified for Neuron by FE GANNIER
	francois.gannier@univ-tours.fr (University of TOURS)
ENDCOMMENT
INCLUDE "Unit.inc"
INCLUDE "Volume.inc"
NEURON {
	SUFFIX inak
	USEION k READ ko WRITE ik
	USEION na READ nai WRITE ina

	RANGE imax, ik, ina, ip
}
PARAMETER {
	Kmna = 40	(mM)
	Kmk = 1		(mM)
	imax = 125	(nA)
}

ASSIGNED {
	v (mV)
	celsius (degC) : 37
	ik (mA/cm2)
	ina (mA/cm2)
	ip (mA/cm2)
	ko (mM)
	nai (mM)
}

BREAKPOINT { 
	ip = (1e-06)*imax/S * (ko/(Kmk + ko))*(nai/(Kmna+nai))

	ina = 3*ip
	ik = -2*ip
} 
