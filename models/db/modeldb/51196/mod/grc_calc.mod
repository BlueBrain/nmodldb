TITLE Cerebellum Granule Cell Model, calcium kinetics
COMMENT  
Reference: E.D'Angelo, T.Nieus, A. Maffei, S. Armano, P. Rossi,
V. Taglietti, A. Fontana, G. Naldi "Theta-frequency bursting and 
resonance in cerebellar granule cells: experimental evidence and 
modeling of a slow K+-dependent mechanism", J. neurosci., 2001,
21,P. 759-770.
ENDCOMMENT

NEURON {
        SUFFIX Calc
        USEION ca READ ica, cao WRITE cai
        RANGE d, beta, cai0
}

UNITS {
        (mV)    = (millivolt)
        (mA)    = (milliamp)
	  (um)    = (micron)
	  (molar) = (1/liter)
        (mM)    = (millimolar)
   	   F      = (faraday) (coulomb)
}

PARAMETER {
        ica             (mA/cm2)
        d = .2          (um)        
        cai0 = 1e-4     (mM)         
        beta = 1.5        (/ms)
}
ASSIGNED{
cao (mM)
}
STATE {
	cai (mM) <1e-04>
}

INITIAL {
        cai = cai0 
}

BREAKPOINT {
       SOLVE conc METHOD derivimplicit
}

DERIVATIVE conc {    
	cai' = -ica/(2*F*d)*(1e4) - beta*(cai-cai0)
}


