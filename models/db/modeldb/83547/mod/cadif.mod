COMMENT
Longitudinal diffusion of sodium (no buffering)
(equivalent modified euler with standard method and
equivalent to diagonalized linear solver with CVODE )
LONGITUDINAL_DIFFUSION DCa*PI*diam*diam/4 {cai}
	 ~ cai << (-fCa*ica/(FARADAY)*PI*diam*(1e4))
ENDCOMMENT

NEURON {
	SUFFIX cadifus
	USEION ca READ cai, ica WRITE cai
	RANGE DCa, cainit, fCa
}

UNITS {
	(molar) = (1/liter)
	(mM) =  (millimolar)
	(um) = (micron)
	(mA) = (milliamp)
	FARADAY = (faraday) (coulomb) 
	PI = (pi) (1)
}

PARAMETER {
	DCa = .6 (um2/ms)
        cainit = 0.000250 (mM)
        fCa = 0.005  (1)
}

ASSIGNED {
	diam  (um)
	ica   (mA/cm2)
		
}

STATE {
	cai (mM) <1e-10>
}

BREAKPOINT {
	SOLVE state METHOD sparse
}

INITIAL{
        cai=cainit
}

KINETIC state {
	COMPARTMENT PI*diam*diam/4 {cai}
	 ~ cai << (-fCa*ica*PI*diam*(1e4)/(2*FARADAY))
}
