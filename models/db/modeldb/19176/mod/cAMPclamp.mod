COMMENT
Mechanism to clamp intracellular cAMP concentration.
ENDCOMMENT

NEURON {
	SUFFIX cAMPclamp
	RANGE del, dur, conc1, conc2
	USEION a  WRITE ai  VALENCE 0
}
UNITS {
        (nA) = (nanoamp)
        (molar) = (1/litre)
}

PARAMETER {
	del (ms)
	dur (ms)      <0,1e9>
        conc1 (molar)   
        conc2 (molar)  
}
ASSIGNED { ai (molar) }

INITIAL {
	ai = conc1
}

BREAKPOINT {
	at_time(del)
	at_time(del+dur)

	if (t < del + dur && t > del) {
		ai = conc2
	}else{
		ai = conc1
	}
}