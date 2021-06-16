TITLE cAMP influx based on cai
COMMENT

	Ca++ binding on second messenger
		p0 (inactive) + nca Ca <-> p1 (active)	; rate cst k1,k2

ENDCOMMENT

UNITS {
    (molar) = (1/liter)
    (mM) = (millimolar)
}

NEURON {
	THREADSAFE
    SUFFIX CN
    USEION ca READ cai
	USEION cn WRITE cni VALENCE 1
    GLOBAL tau
}

PARAMETER {
	: all values can be adjusted in hoc files
	tau_ca = 1500	(ms)
	tau = 6000		(ms)
	cn_init = 2e-5	(mM)
	kD = 0.05		(1)
	minca = 1.5e-4	(mM)
}

ASSIGNED {
    cai (mM)
    cinf (mM)
}

STATE {
	cni (mM)
}

INITIAL {
    cni = cn_init
}


BREAKPOINT {
	SOLVE state METHOD derivimplicit
}


DERIVATIVE state { LOCAL dcn

	influx(cai)
	:dcn = cinf/tau_ca - cni/tau
	cni' = cinf/tau_ca - cni/tau
}


PROCEDURE influx( cai(mM) ) {
	:TABLE cinf DEPEND minca, kD
	:      FROM 0 TO 0.02 WITH 1000
    
	if (cai < minca) {
		cinf = 0
	} else {
		cinf = (cai-minca)*kD
	}
}





