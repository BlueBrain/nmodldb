TITLE Calcium ion accumulation without diffusion and buffering
 
                 
NEURON {
	SUFFIX cabalan
	USEION ca READ cai, ica WRITE cai
	RANGE  cainit, fCa , icapump,icapumpmax,km
}

UNITS {
	(molar) = (1/liter)
	(mM) =  (millimolar)
	(um) =  (micron)
	(mA) =  (milliamp)
	FARADAY = (faraday) (coulomb) 
	PI = (pi) (1)
}

PARAMETER {
         fCa = 0.05  (1)
         cainit = 0.00002 (mM)
        dt    (ms)
        celsius = 35  (degC)
        icapumpmax  = 0.00191  (mA/cm2)
        km = 0.000500         (mM)
         }

ASSIGNED {
	diam  (um)
	ica   (mA/cm2)
          icapump (mA/cm2)
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
        icapump = icapumpmax*(1/(1 + km/cai))
	COMPARTMENT PI*diam*diam/4 {cai}
	 ~ cai << (-fCa*(ica +icapump)*PI*diam*(1e4)/(2*FARADAY))
}
