TITLE Potassium ion accumulation wo cleft
COMMENT
	modified From DiFrancesco & Noble 1985 Phil Trans R Soc Lond 307:353-398 
    modified for Neuron by FE GANNIER
	francois.gannier@univ-tours.fr (University of TOURS)
ENDCOMMENT
INCLUDE "Unit.inc"
INCLUDE "Volume.inc"
NEURON {
	SUFFIX K_acc
	USEION k READ ki, ik, ko WRITE ki, ko
	RANGE ik, pf, kb
}

PARAMETER {
	ik 			(mA/cm2)
	kb = 4		(mM)
	
	pf = 0.7	(/s)
}

STATE {
	ki START 140	(mM)
	ko START 4		(mM)
}

LOCAL ViF, VeF
INITIAL {
	VERBATIM
		ki = _ion_ki;
		ko = _ion_ko;
	ENDVERBATIM
	VeF = (1e-3)*Ve*F/S
	ViF = (1e-3)*Vi*F/S
}

BREAKPOINT {
:	SOLVE state METHOD derivimplicit
	SOLVE state METHOD cnexp
}

: Rq  Ko = Kc
DERIVATIVE state { 
	
	ki' = -ik /ViF
	ko' = ik /VeF -(0.001)*pf*(ko - kb)
}
