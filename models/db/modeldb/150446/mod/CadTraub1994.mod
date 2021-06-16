TITLE Calcium dynamics in a submembrane shell

COMMENT
  from paper "A branching dendritic model of a rodent CA3 pyramidal neurone." Traub RD et al. J Physiol. (1994) 
  implemented by Nikita Vladimirov <nikita.vladimirov@gmail.com>
ENDCOMMENT

NEURON {
	SUFFIX Cad
	USEION ca READ ica WRITE cai
	RANGE  phi, beta
}

UNITS {
	(mA)	= (milliamp)
 	(molar) = (1/liter)			: moles do not appear in units
 	(mM)	= (millimolar)
}

PARAMETER {
	phi	 	
	beta 	(/ms)
}

STATE {	cai (mM) }

INITIAL { 
	cai = 0
}

ASSIGNED { 
	ica		(mA/cm2) 
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	if( cai < 0 ){ cai = 0 }
}

DERIVATIVE state { 
	cai' = - phi * ica - beta * (cai)
}
