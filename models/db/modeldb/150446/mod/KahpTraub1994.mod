TITLE Potassium after-hyperpolarization, K(AHP) current

COMMENT
  from Table 3 of paper "A branching dendritic model of a rodent CA3 pyramidal neurone." Traub RD et al. J Physiol. (1994) 
  implemented by Nikita Vladimirov <nikita.vladimirov@gmail.com>
ENDCOMMENT

NEURON { 
	SUFFIX Kahp
	USEION k READ ek WRITE ik
	USEION ca READ cai
	RANGE gbar, g, i
}

UNITS { 
	(mV) = (millivolt) 
	(mA) = (milliamp) 
 	(molar) = (1/liter)			: moles do not appear in units
 	(mM)	= (millimolar)
} 

PARAMETER { 
	gbar = 0.0 	(S/cm2)
}
 
ASSIGNED { 
	v		(mV)
	cai		(mM)
	ek		(mV)
	ik 		(mA/cm2) 
	i 		(mA/cm2) 
	g		(S/cm2)
	minf
	mtau    (ms) 
}
 
STATE {	m }

BREAKPOINT { 
	SOLVE states METHOD cnexp
	g = gbar * m 
	i = g * ( v - ek ) 
	ik = i
}
 
INITIAL { 
	rates( v, cai )
	m = minf
}
 
DERIVATIVE states { 
	rates( v, cai )
	m' = (minf - m) / mtau
}

PROCEDURE rates( v (mV), ca (mM) ) { 
	LOCAL  a, alpham, betam
	UNITSOFF 
	a = 0.2 * 1e-4 * ca
	if( a < 0.01 ) {
		alpham = a
	} else {
		alpham = 0.01
	}
	betam = 0.001
	minf   = alpham / ( alpham + betam )
	mtau   = 1 / ( alpham + betam )
	UNITSON
}
