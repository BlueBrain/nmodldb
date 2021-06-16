TITLE intracellular Calcium concentration dynamic
: from BEELER & REUTER, J.Physiol, 1977

COMMENT
	 Cardiac intracellular calcium 
ENDCOMMENT

NEURON {
	THREADSAFE
	SUFFIX Cadynam
	USEION ca READ ica, cai WRITE cai
	USEION cs READ ics VALENCE 2
	RANGE ics
}

UNITS {
	(mM) = (milli/liter)
	(mA) = (milliamp)
}

ASSIGNED{
	ica		(mA/cm2)
	ics		(mA/cm2)
}

PARAMETER {
}

STATE {
	cai START 0.001 (mM)  <1e-6>
}

INITIAL {
	VERBATIM
		cai = _ion_cai;
	ENDVERBATIM
}

BREAKPOINT {
	SOLVE state METHOD derivimplicit
}

DERIVATIVE state { 
	cai' = -1e-7*ics+0.07*(1e-7-cai)
}
