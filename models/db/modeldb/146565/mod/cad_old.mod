

NEURON {
	SUFFIX cad_old
	USEION ca READ ica, cai WRITE cai
	GLOBAL taur,cainf
}

CONSTANT {
	FARADAY = 96489				
}

PARAMETER {
	taur	= 100				
	cainf	= 0.0001
}
ASSIGNED {
	ica
}
STATE {
	cai	 
}

INITIAL {
	cai = cainf
}

BREAKPOINT {
	SOLVE state METHOD euler
}

DERIVATIVE state {
	cai'=-(100000)*ica/(2 * FARADAY )-(cai-cainf)/taur
}
