: $Id: kahpRT03.mod,v 1.1.1.1 2005/12/15 15:16:39 hines Exp $
TITLE Potasium AHP type current for RD Traub, J Neurophysiol 89:909-921, 2003

COMMENT

	Implemented by Maciej Lazarewicz 2003 (mlazarew@seas.upenn.edu)

ENDCOMMENT

UNITS { 
	(mV) = (millivolt) 
	(mA) = (milliamp) 
} 

NEURON { 
	SUFFIX kahpRT03
	USEION k READ ek WRITE ik
	USEION ca READ cai
	RANGE gmax, ik, g, i
        GLOBAL scale
}

PARAMETER { 
	gmax = 0.0 	(mho/cm2)
	v		(mV) 
	ek 		(mV)  
	cai		(1)
        scale = 1e-3
}
 
ASSIGNED { 
  cas
  g 
  i
  ik 		(mA/cm2) 
  alpha beta	(/ms)
}
 
STATE {
	m
}

BREAKPOINT { 
	SOLVE states METHOD cnexp
        g = gmax * m
	i = g * ( v - ek ) 
        ik = i
}
 
INITIAL { 
	rates()
	m = alpha / ( alpha + beta )
	m = 0
}
 
DERIVATIVE states { 
	rates()
	m' = alpha * ( 1 - m ) - beta * m 
}

UNITSOFF 

PROCEDURE rates() { 

  cas=cai*scale
  if (cas < 100) {
    alpha = cas / 10000
  }else{
    alpha = 0.01
  }
  beta = 0.01
}

UNITSON
