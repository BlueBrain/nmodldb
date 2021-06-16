TITLE Time-independent (delayed) K+ current Ik1
COMMENT
	modified From DiFrancesco & Noble 1985 Phil Trans R Soc Lond 307:353-398 
    modified for Neuron by FE GANNIER
	francois.gannier@univ-tours.fr (University of TOURS)
ENDCOMMENT
INCLUDE "Unit.inc"
INCLUDE "Volume.inc"
NEURON {
	SUFFIX ik1
	USEION k READ ek, ko WRITE ik
	RANGE g, ik
}

PARAMETER {
	g = 920		(uS)
	Km1 = 210	(mM)
}

ASSIGNED {
	v (mV)
	celsius (degC) : 37
	ik (mA/cm2)
	minf 
	mtau (ms)  
	ek (mV)
	ko (mM)
	ki (mM)
}

LOCAL RT
INITIAL {
	RT = (1000)*R*(273.15+celsius)
}

BREAKPOINT { 
	ik = (1e-06)* g/S * (ko/(ko + Km1))*((v-ek)/(1 + exp((v - ek + 10)*2*F/RT)))
}
