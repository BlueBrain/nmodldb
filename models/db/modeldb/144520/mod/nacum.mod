TITLE Sodium ion accumulation
COMMENT
	modified From DiFrancesco & Noble 1985 Phil Trans R Soc Lond 307:353-398 
    modified for Neuron by FE GANNIER
	francois.gannier@univ-tours.fr (University of TOURS)
ENDCOMMENT
INCLUDE "Unit.inc"
INCLUDE "Volume.inc"
NEURON {
	SUFFIX Na_acc
	USEION na READ ina, nai  WRITE nai
	RANGE ina
}

PARAMETER {
	ina			(mA/cm2)
}

STATE {
	nai START 8	(mM)
}

LOCAL ViF
INITIAL {
	VERBATIM
		nai = _ion_nai;
	ENDVERBATIM
	ViF = (1e-3)*Vi*F/S
}

BREAKPOINT {
	SOLVE state METHOD derivimplicit
}

DERIVATIVE state {
	nai' = -ina/(ViF)
}
