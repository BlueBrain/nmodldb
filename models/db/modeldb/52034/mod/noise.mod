
COMMENT
current noise injection, changes between + and - imax every dt
Paul Bush 1995
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS noise
	NONSPECIFIC_CURRENT i
	RANGE imax
:	GLOBAL seed
}

ASSIGNED {
	rn
}

UNITS {
	(nA) = (nanoamp)
}

PARAMETER {
	imax=1	        (umho)
:	seed=1 
}

INITIAL {
:	srandom(seed)
	rn = (1/(2^31))*2
}

ASSIGNED { i (nA) }

BREAKPOINT {

SOLVE dum	:has to be in a proc otherwise there is an error

}

PROCEDURE dum() {
i = (scop_random()*rn-1)*imax
}

