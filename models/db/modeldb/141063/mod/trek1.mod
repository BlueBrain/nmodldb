TITLE voltage- and temperature-sensitive Trek-1 channel
: a.k.a. KCNK2, K2p2.1
: *** voltage-independent mode; only temperature dependent ***
: Only steady state, no kinetics at all.
:
: written by Patricio Orio to resemble the figures from: 
: Bockenhauer et al. (Nature Neurosci. 4(5):486-491(2001)), 
: Kang et al. (J Physiol 564(1):103-116(2005)) and 
: Maingret et al. (EMBO J. 19(11):2483-2491(2000)) 

NEURON {
	SUFFIX trek1
	USEION k READ ek WRITE ik
	RANGE gkbar
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gkbar   = .005  (mho/cm2)
	Pmin = 0.1  (1)    :Temperature-independent leak
	q10 = 5   (1)
	Tm  =  39   (degC)
}   

ASSIGNED {
	ek      (mV)
	celsius	(degC)
	ik      (mA/cm2)
	v	(mV)
	o
}

BREAKPOINT {
	Ps()
    ik  = o * gkbar * (v - ek)
}

INITIAL {
    o = Pmin + (1-Pmin)/(1+exp(-(celsius-Tm)*q10/10 (degC))) 
	ik  = o * gkbar * (v - ek)
}

PROCEDURE Ps(){
    o = Pmin + (1-Pmin)/(1+exp(-(celsius - Tm)*q10/10 (degC))) 
}