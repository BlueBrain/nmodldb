COMMENT
I added area in LON.._DIF.. expression and buffering

Longitudinal diffusion of sodium (with buffering coefficient f)
(equivalent modified euler with standard method and
equivalent to diagonalized linear solver with CVODE )
	LONGITUDINAL_DIFFUSION D*PI*diam*diam/4 {nai}
ENDCOMMENT

NEURON {
	SUFFIX nadifl
	USEION na READ ina WRITE nai
	RANGE D,nainit,f
}

UNITS {
	(mM) = (milli/liter)
	(um) = (micron)
	FARADAY = (faraday) (coulomb)
	PI = (pi) (1)
}

PARAMETER {
	D = .6 (um2/ms)
        nainit = 3.54 (mM) 
        f = 1.00
}

ASSIGNED {
	ina (milliamp/cm2)
	diam (um)
}

STATE {
	nai (mM) <1e-4>
}

BREAKPOINT {
	SOLVE state METHOD sparse
}

INITIAL{
        nai=nainit
}

KINETIC state {
	COMPARTMENT PI*diam*diam/4 {nai}
	~ nai << (-f*ina*PI*diam*(1e4)/(FARADAY))
}

