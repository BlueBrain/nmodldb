COMMENT

I P/Q from Bischofberger et al. J neuroscience 2002

ENDCOMMENT

NEURON {
	SUFFIX capq
	USEION ca READ eca WRITE ica
	RANGE gca, gbar, i
	RANGE minf, mtau : would be OK for these to be GLOBAL
	GLOBAL midv, mslope, k1, k2, v1, v2, temp
}


UNITS { : units that are not in the units database should be declared here
  (mV) = (millivolt)
  (mA) = (milliamp)
  (uA) = (microamp)
  (S) = (siemens)
}

PARAMETER {
	: set to the values described in the aforementioned paper
	eca = 110 (mV) : this value will have no effect. set in hoc code
	gbar = 1 (S/cm2) 	
	midv = -3.9 (mV)
	mslope = 7.1 (mV)
	k1= 0.22 : this value has been modified from 1.12 to 0.22 in order to observe both d-ADF and h-ADF in the model
	k2= 0.14
	v1= 31.5 (mV)
	v2= 8.6 (mV)
	delay=0.2 (ms)
	temp=0
}

ASSIGNED {
	: either assigned by the system (e.g., v and i) or by us
	v (mV)
	i 	(mA/cm2)
	ica 	(mA/cm2)
	gca	(S/cm2)
	minf
	mtau (ms)
}

STATE { m }

INITIAL { 
	rates(v)
	m = minf
}

BREAKPOINT {
		SOLVE states METHOD cnexp
        	gca = gbar*m
		i = gca*(v - eca) : for convenience, "i" is declared as range so that it can be studied as a seperate current coming from this mechanism.
		ica = i
}

DERIVATIVE states {  
		rates(v)
        m' = (minf - m)/mtau
}

PROCEDURE rates(v (mV)) 
{
	minf = 1/ ( 1 + exp( (midv-v)/mslope ) )
	mtau = ((2.3^-temp)/(k1*exp(v/v1)+k2*exp(-v/v2)))+delay
}


