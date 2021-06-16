TITLE Cerebellum Granule Cell Model, Leakage

COMMENT   
Reference: E.D'Angelo, T.Nieus, A. Maffei, S. Armano, P. Rossi,
V. Taglietti, A. Fontana, G. Naldi "Theta-frequency bursting and 
resonance in cerebellar granule cells: experimental evidence and 
modeling of a slow K+-dependent mechanism", J. neurosci., 2001,
21,P. 759-770.
ENDCOMMENT
 
NEURON { 
	SUFFIX GrC_Lkg1 
	NONSPECIFIC_CURRENT il
	RANGE el, gl,i
} 
 
UNITS { 
	(mA) = (milliamp) 
	(mV) = (millivolt) 
} 
 
PARAMETER { 
		gl = 5.68e-5 (mho/cm2)
	      el = -58 (mV)
} 

ASSIGNED {
      v (mV) 
      il (mA/cm2) 
	i (mA/cm2) 
}
  
BREAKPOINT { 
	il = gl*(v - el) 
	i = il
} 

