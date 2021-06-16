TITLE Sodium ion accumulation
: Intracellular sodium ion accumulation 

NEURON {
	THREADSAFE
	SUFFIX Na_acc
	USEION na READ ina,nai WRITE nai
	RANGE Vi, Naneutral : electroneutral sodium accumulation
}

UNITS {
	(mV) = (millivolt)
	(um) = (micron)
	(mM) = (milli/liter)
	(mA) = (milliamp)
	F = (faraday) (coulombs)
}

PARAMETER {
	Naneutral = 3e-5   (mA/cm2)   <0,1e6>
	Vi = 1.3668e-08     (cm3)
	ina		   (mA/cm2)
}

STATE {
	nai START 9.19	(mM)
}

LOCAL ViF
INITIAL {
	VERBATIM
	nai = _ion_nai;
	
	ENDVERBATIM
	ViF = F*Vi*2e4
}

BREAKPOINT {
	SOLVE state METHOD derivimplicit
}

DERIVATIVE state {
	nai' = -(ina-Naneutral)/(ViF) 
}

