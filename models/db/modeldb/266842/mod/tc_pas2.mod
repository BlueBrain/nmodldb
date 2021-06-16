TITLE Thalamocortical neuron passive leakage channels

COMMENT

	Leakage Na+ and K+ currents responsible for resting potential
	Implementation of Meijer et al., 2011
	Written by Xu Zhang, UConn, 2018
	
ENDCOMMENT

NEURON {
	SUFFIX tcpas2
	USEION na READ ena WRITE ina
	USEION k READ ek WRITE ik
	NONSPECIFIC_CURRENT i
	RANGE g_kl, g_nl
}


UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
	(nA) = (nanoamp)
	(pA) = (picoamp)
	(S)  = (siemens)	
}

PARAMETER {
	g_kl  = 5e-5  (mho/cm2)
	g_nl  = 2e-5  (mho/cm2)
	ena     = 45    (mV)
	ek      = -95   (mV)
	v               (mV)
}

ASSIGNED {
	i (mA/cm2)
 	ik (mA/cm2)
	ina (mA/cm2)
}

BREAKPOINT {
	ik = g_kl * (v - ek)
	ina = g_nl * (v - ena)
	i = ik + ina
}