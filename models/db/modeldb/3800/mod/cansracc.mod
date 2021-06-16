COMMENT
	calcium accumulation into NSR, taken from from Courtemanche et al Am J Physiol 1998 275:H301
ENDCOMMENT

NEURON {
	SUFFIX ca_nsracc
	USEION cu READ cui WRITE cui VALENCE 2 : calcium ion in nsr
	USEION cr READ  cri  VALENCE 2 : calcium ion in jsr
	USEION ca READ cai
	RANGE Vrel, Vup, Iupmax, Kup, cupmax
}

UNITS {
	(mM) = (milli/liter)
	(mA) = (milliamp)
	
}
ASSIGNED{
	cri 			(mM)
	cai 			(mM)
}
PARAMETER {
	Vrel =  96.48e-12     	(cm3)
	Vup = 1109.52e-12  	(cm3) 
	Kup = 0.00092          	(mM) 
	Iupmax = 0.005   	(mM/ms)
	cupmax = 15		(Mm)
	
	
}

STATE {
	cui START 1.5	(mM)
}

INITIAL {
	VERBATIM
	cui = _ion_cui;
	
	ENDVERBATIM
}

BREAKPOINT {
	SOLVE state METHOD derivimplicit
}

DERIVATIVE state {
	cui' = Iupmax/(1 + (Kup/cai)) - cui*Iupmax/cupmax - (cui - cri)/180*Vrel/Vup
}
