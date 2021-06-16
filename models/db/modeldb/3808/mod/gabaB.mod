TITLE GABAb synapse

COMMENT
        simple alpha-synapse that generates a single PSP   
        *********************************************
        reference:   	McCormick, Wang & Huguenard (1993) 
			Cerebral Cortex 3(5), 387-398
        found in:       cat reticular nucleus of thalamus
        *********************************************
	Assembled for MyFirstNEURON by Arthur Houweling
ENDCOMMENT
					       
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS GABAbSynapse
	USEION k READ ek 
	RANGE onset, gmaxIPSP, e, g, i, w
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) 	= (nanoamp)
	(mV)	= (millivolt)
	(nS) 	= (nanomho)
}

PARAMETER {
	onset= 25	(ms)
	gmaxIPSP= 0	(nS)
	w= 1				: weight factor for gmaxIPSP
	kf= 1				: kinetic factor
	ek		(mV)
	v		(mV)
	celsius		(degC)
}

ASSIGNED { 
	i 		(nA)  
	g 		(nS)
	tadj
}

UNITSOFF
INITIAL {
        tadj = 3^((celsius-23.5)/10)
}    

BREAKPOINT { LOCAL tt
	tt = (t-onset)*tadj
	if (t>onset) {
	  g = w*gmaxIPSP * (0.84*exp(-tt/283) + 
 	      0.16*exp(-tt/1026)) * ((1-exp(-tt/112))^4)/0.31
	}
	else {g = 0}
	i = g * (v-ek)
}
UNITSON

