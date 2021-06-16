COMMENT
	cardiac calcium accumulation into jsr, from Courtemanche et al Am J Physiol 1998 275:H301
ENDCOMMENT

NEURON {
	SUFFIX ca_jsracc
	USEION ca READ cai
	USEION cu READ cui VALENCE 2 : calcium ion in nsr
	USEION cr READ icr, cri WRITE cri VALENCE 2 : calcium ion in jsr
	RANGE Csqn, Kmcsqn, Iupmax, Kup
}

UNITS {
	(mM) = (milli/liter)
	(mA) = (milliamp)
	
}
ASSIGNED{
	icr 		(mA/cm2)
	cui		(mM)
	cai		(mM)
}
PARAMETER {
	Csqn =  10      (mM)
	Kmcsqn = 0.8 	(mM)
	Iupmax = 0.005   (mM/ms)
	Kup = 0.00092    (mM)  
	
	
}
 

STATE {
	cri START 1.5	(mM)
}
INITIAL {
	VERBATIM
	cri = _ion_cri;
	
	ENDVERBATIM
}

BREAKPOINT {
	SOLVE state METHOD derivimplicit
}

DERIVATIVE state {
	cri' = ((cui - cri)/180 - 0.5*icr)/(1 + Csqn*Kmcsqn/((cri + Kmcsqn)*(cri + Kmcsqn)))
}
