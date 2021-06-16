NEURON { SUFFIX ibackg : cardiac (atrial) background currents from Courtemanche et al Am J Physiol 1998 275:H301
	USEION na READ ena WRITE ina
	USEION ca READ eca WRITE ica
	RANGE gbNa, gbCa
	GLOBAL dummy : prevent vectorization for use with CVODE
}

UNITS { (mV) = (millivolt)  
	(mA) = (milliamp) 
	(mM) = (milli/liter) 
	
}

PARAMETER {
	gbNa =  0.001348e-3  (S/cm2)  <0,1e9>
	gbCa =  0.00226e-3  (S/cm2)    <0,1e9>
}

ASSIGNED {
	celsius (degC)
	v (mV) 
	ina (mA/cm2)  
	ica (mA/cm2)
	ena (mV)
	eca (mV)
	dummy
}

BREAKPOINT {
	
	ina = gbNa*(v - ena)

	ica = gbCa*(v - eca)
}

