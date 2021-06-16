TITLE sodium potassium pump
:  Adapted from Leo Medina's implementation from Lindblad et al Am J Physiol 1996 275:H1666

: Original model has been modified to assume constant nai

NEURON {
	SUFFIX aNaKpump
	USEION k READ ko WRITE ik	
	USEION na READ nai WRITE ina
	RANGE INaKmax, ina, ink, Kmko, Kmnai, ik
	GLOBAL dummy : prevent vectorization for use with CVODE
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
  	(mM) = (millimolar)
	
}

PARAMETER {
	INaKmax = 0.009726135 (mA/cm2) <0,1e6>
	Kmnai = 5.46 (mM)    <0,1e6>
	Kmko = 0.621 (mM)    <0,1e6>
	Q10NaK = 1.16
}

ASSIGNED {
	celsius (degC)
	v (mV)
	ko (mM)
	nai (mM)
	ik (mA/cm2)
	ina (mA/cm2)
	ink (mA/cm2)
	dummy
}

BREAKPOINT { LOCAL fnk
	
	fnk = (v + 150)/(v + 200)
				
	ink = INaKmax*fnk*((nai/(nai+Kmnai))^3)*((ko/(ko+Kmko))^2) : Changed this line to reflect the exponents given in Schild 1994, instead of the orginal exponents in Leo's model.
	
	if (celsius >= 37) {
		ink=Q10NaK*ink
	}
	
	ina = 3*ink
	ik = -2*ink
}