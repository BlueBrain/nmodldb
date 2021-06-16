TITLE Ih current from HCN1 units

COMMENT

NEURON implementation of a HCN-channel
Kinetical Scheme: Hodgkin-Huxley (n)

Modified from Khaliq et al. J. Neurosci. 23(2003)4899
Model HCN includes the calculation of the gating current

Reference: Akemann et al., Biophys. J. (2009) 96: 3959-3976

Laboratory for Neuronal Circuit Dynamics
RIKEN Brain Science Institute, Wako City, Japan
http://www.neurodynamics.brain.riken.jp

Date of Implementation: April 2007
Contact: akemann@brain.riken.jp

ENDCOMMENT

NEURON {
	SUFFIX HCN
	NONSPECIFIC_CURRENT i
	RANGE i, igate, gbar, eh
	GLOBAL ninf, taun
	GLOBAL gateCurrent, gunit
}

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
	(nA) = (nanoamp)
	(pA) = (picoamp)
	(S)  = (siemens)
	(nS) = (nanosiemens)
	(pS) = (picosiemens)
	(um) = (micron)
	(molar) = (1/liter)
	(mM) = (millimolar)	
}

CONSTANT {
	e0 = 1.60217646e-19 (coulombs)
	q10 = 2.7
	
	cvn = 90.1 (mV)
	ckn = -9.9 (mV)
		
	cct = 0.19 (s)
	cat = 0.72 (s)
	cvt = 81.5 (mV)
	ckt = 11.9 (mV)

	zn = 2.5691 (1)		: valence of n-gate
}

PARAMETER {
	gateCurrent = 0		: gating currents ON = 1 OFF = 0

	gbar = 0.0001 (S/cm2)
	gunit = 0.68 (pS)
	eh = -30 (mV)
}

ASSIGNED {
	celsius (degC)
	v (mV)
	
	i (mA/cm2)
	igate (mA/cm2)
	qt (1)
	nc (1/cm2)

	ninf (1)
	taun (ms)
}

STATE { n }

INITIAL {
	nc = (1e12) * gbar / gunit
	qt = q10^((celsius-22 (degC))/10 (degC))
	rates(v)
	n = ninf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	i = gbar * n * (v - eh)
	igate = - nc * (1e6) * e0 * zn * ngateFlip()

	if (gateCurrent != 0) { 
		i = i + igate
	}
}

DERIVATIVE states {
	rates(v)
	n' = (ninf-n)/taun
}

PROCEDURE rates(v (mV)) {
	ninf = 1 / ( 1+exp(-(v+cvn)/ckn) )
	taun = (1e3) * ( cct + cat * exp(-((v+cvt)/ckt)^2) ) / qt
}

FUNCTION ngateFlip() (1/ms) {
	ngateFlip = (ninf-n)/taun
}


