COMMENT
	Internal Ca concentration from Av-Ron and Vidal 1999
	Implemented by C. Weaver, 2003
ENDCOMMENT

UNITS {
	(molar) =	(1/liter)
	(mM) =	(millimolar)
	(mA) = (milliamp)
	(mV) = (millivolt)
	FARADAY = (faraday) (kilocoulombs)
	R	= (k-mole) (joule/degC)

}

NEURON {
	SUFFIX cad
	USEION ca READ ica, cai, eca WRITE cai
	GLOBAL Kp, Rca, cainf
    RANGE cai, cao
}

PARAMETER {
	celsius 	(degC)
        v       (mV)
	: Rca=5
	: Kp=1
	Kp = 0.0005	(mM-cm2/mA-ms)
	Rca = 0.003	(/ms)
	cainf	= 100e-6(mM)
	: cainf	= 100e-5(mM)
	: cai	(mM)
}


STATE {
	cai (mM)
}

INITIAL {
	cai = cainf
}

ASSIGNED {
	ica 	(mA/cm2)
      ena      (mV)
	eca	(mV)
}

BREAKPOINT {
	SOLVE state METHOD cnexp
}

DERIVATIVE state {     : exact when v held constant; integrates over dt step
	cai' = Kp * (-1 * ica) - Rca*cai
}

