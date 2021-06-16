COMMENT

	Maciej T. Lazarewicz, mlazarew@seas.upenn.edu

ENDCOMMENT



TITLE rcadecay

INDEPENDENT {t FROM 0 TO 1 WITH 10 (ms)}

NEURON {

	SUFFIX cadecay
	USEION ca READ ica WRITE cai
	RANGE  phi, beta
}

UNITS {

	(molar) = (1/liter)
	(mM)	= (millimolar)
	(um)	= (micron)
	(mA)	= (milliamp)
}

PARAMETER {

	phi             = 0.13e3
	beta            = 0.075
}
ASSIGNED {

  ica (milliamp/cm2)
}

STATE {	cai (milli/liter) }

INITIAL {
 
  cai= - phi * ica/ beta
  :cai = 0.2
}

BREAKPOINT {

  if       ( cai < 0 )      { cai = 0 
  } else                    { SOLVE state METHOD cnexp }

}

UNITSOFF

DERIVATIVE state { cai' = - phi * ica - beta * cai  }

UNITSON