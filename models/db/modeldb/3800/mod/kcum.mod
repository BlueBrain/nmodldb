TITLE Potassium ion accumulation
: Intracellular potassium ion accumulation 

NEURON {
	SUFFIX K_acc
	USEION k READ ki, ik WRITE ki
	RANGE Vi, Kneutral: electroneutral K accumulation
}

UNITS {
	
	(mV) = (millivolt)
	(mM) = (milli/liter)
	(mA) = (milliamp)
	F = (faraday) (coulombs)
}

PARAMETER {

	Vi = 13668e-12 (cm3)
	ik 	(mA/cm2)
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
	
