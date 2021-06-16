TITLE Cai.mod    

COMMENT
Intracellular Ca pool

Author: Fredrik Edin, 2003
Address: freedin@nada.kth.se

ENDCOMMENT
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
	(mM) = (milli/liter)
}

NEURON {
        SUFFIX Cai
        USEION ca READ ica WRITE cai
        RANGE Alpha, tau
}
 
PARAMETER {
	Alpha 	= 0.000667	(mM-cm2/ms-mA) 
	tau 	= 240 		(ms)
}

ASSIGNED {
	gca (mho/cm2)
        ica (mA/cm2)
}

STATE {
	cai (mM) <1e-12>
}
 
BREAKPOINT { SOLVE state METHOD cnexp }

DERIVATIVE state {
	cai' = - Alpha * ica - cai / tau
}