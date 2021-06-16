TITLE Potassium ion accumulation
: Intracellular potassium ion accumulation 

NEURON {
	THREADSAFE
	SUFFIX K_acc
	USEION k READ ki, ik WRITE ki
	RANGE Vi, ik, Kneutral: electroneutral K accumulation
}

UNITS {
	
	(mV) = (millivolt)
	(mM) = (milli/liter)
	(mA) = (milliamp)
	F = (faraday) (coulombs)
}

PARAMETER {
	Vi = 1.3668e-08 (cm3)
	ik = 0.000445123 (mA/cm2)
	Kneutral=3e-5 (mA/cm2) <0,1e6>
}

STATE {
	ki START 141.13	(mM)
}

LOCAL ViF
INITIAL {
	VERBATIM
		ki = _ion_ki;
	ENDVERBATIM
	ViF = Vi*F*2e4
}

BREAKPOINT {
	SOLVE state METHOD derivimplicit
}

DERIVATIVE state {
	ki' = -(ik-Kneutral) /(ViF)
}
