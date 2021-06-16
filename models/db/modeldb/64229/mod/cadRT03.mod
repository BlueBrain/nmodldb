: $Id: cadRT03.mod,v 1.1.1.1 2005/12/15 15:16:39 hines Exp $
TITLE Calcium dynamics for RD Traub, J Neurophysiol 89:909-921, 2003

COMMENT

	Implemented by Maciej Lazarewicz 2003 (mlazarew@seas.upenn.edu)

ENDCOMMENT

NEURON {
	SUFFIX cadRT03
	USEION ca READ ica WRITE cai
	RANGE  phi, beta
	GLOBAL ceiling
}

UNITS {
	(mA)	= (milliamp)
}

PARAMETER {
	phi = 26000	(1)
	beta = .02	(/ms)
	ceiling	= 1000	(1)
}

STATE {	cai (1) }

INITIAL { 
	cai = 0.0 
}

ASSIGNED { 
	ica		(mA/cm2) 
}
	
BREAKPOINT {
  SOLVE state METHOD cnexp
  if( cai<0 ){ printf("cai less than 0 at %g (%g)\n",t,ica) }
  ica=0
}

DERIVATIVE state { 
	cai' = - phi * ica - beta * cai 
}
