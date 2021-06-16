COMMENT
Slavish transliteration of equation from Purvis & Butera 2005.
Not a remapping of their Ca accum mechanism in a form that would use density parameters.

Cannot be used outside of the context of their model 
(K1 implicitly assumes a particular volume and surface/volume ratio, which 
depend on cell geometry).

In the absence of calcium currents, cai will monoexponentially decay toward 0.
ENDCOMMENT

NEURON {
	SUFFIX cax
	USEION ca READ ica WRITE cai
	RANGE K1, K2
}

UNITS {
	(mM) = (milli/liter)
	(mA) = (milliamp)
	F = (faraday) (coulombs)
}

PARAMETER {
	K1 = 0.02 (mM-cm2/ms-mA)
	K2 = 0.04 (/ms)
	cai0 = 0.0000604 (mM)
}

ASSIGNED {
	 ica (mA/cm2)
}

STATE {
	cai (mM)
}

INITIAL {
	cai = cai0
}

BREAKPOINT {
	SOLVE integrate METHOD derivimplicit
}

DERIVATIVE integrate {
	cai' = -(K1*ica + K2*cai)
}
