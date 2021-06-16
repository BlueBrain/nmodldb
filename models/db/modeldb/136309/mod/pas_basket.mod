TITLE Passive current meant to mimic Basket cells on the soma

COMMENT

	Implemented by Erin Munro
ENDCOMMENT

INDEPENDENT { t FROM 0 TO 1 WITH 1 (ms) }

UNITS { 
	(mV) = (millivolt) 
	(mA) = (milliamp) 
  (S) = (siemens)
} 
NEURON { 
	SUFFIX pasBasket
	NONSPECIFIC_CURRENT i
	RANGE g, e
}
PARAMETER { 
	g = 0.0 	   (mho/cm2)
	e = -70		   (mV)  
} 
ASSIGNED {
    v (mV)
    i (mA/cm2)
}
BREAKPOINT { 
	i = g*(v - e) 
} 
INITIAL {
    :g = 0
    :e = -50
}

