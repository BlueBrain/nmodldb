TITLE Cerebellum Granule Cell Model, Gaba A leakage
COMMENT  
Reference: E.D'Angelo, T.Nieus, A. Maffei, S. Armano, P. Rossi,
V. Taglietti, A. Fontana, G. Naldi "Theta-frequency bursting and 
resonance in cerebellar granule cells: experimental evidence and 
modeling of a slow K+-dependent mechanism", J. neurosci., 2001,
21,P. 759-770.
ENDCOMMENT
 
NEURON { 
	SUFFIX GrC_Lkg2 
	NONSPECIFIC_CURRENT il
	RANGE egaba, ggaba , i
} 
 
UNITS { 
	(mA) = (milliamp) 
	(mV) = (millivolt) 
} 
 
PARAMETER { 
		ggaba = 2.17e-5 (mho/cm2)
	      egaba = -65 (mV)
} 

ASSIGNED { 
       v (mV) 
       il (mA/cm2) 
	 i (mA/cm2) 
}

BREAKPOINT { 
	il = ggaba*(v - egaba) 
	i =il
} 

