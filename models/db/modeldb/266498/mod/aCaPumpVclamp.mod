:Calcium Pump
:  Adapted from Leo Medina's implementation from Lindblad et al Am J Physiol 1996 275:H1666

: Original model has been modified to assume constant nai

NEURON {
	SUFFIX aCaPump
	USEION ca READ cai WRITE ica	
	RANGE ICaPmax, KmCa, ica
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
  	(mM) = (millimolar)
	
}

PARAMETER {
	ICaPmax = 0.000859437(mA/cm2) <0,1e6> 
	KmCa = .0005 (mM)    <0,1e6> 
	Q10CaP = 2.30
}

ASSIGNED {
	celsius (degC)
	v (mV)
	cai (mM)
	ica (mA/cm2)
	icap (mA/cm2)
}

BREAKPOINT {
	
	icap = ICaPmax*(cai/(cai+KmCa)) 
	
	if (celsius >= 37) {
		icap=Q10CaP*icap
	}
	ica=icap
}