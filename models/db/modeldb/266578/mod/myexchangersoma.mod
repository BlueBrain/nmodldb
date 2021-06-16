TITLE sodium calcium exchange
: FORREST MD (2014) Two Compartment Model of the Cerebellar Purkinje Neuron
: taken from Courtemanche et al Am J Physiol 1998 275:H301

NEURON {
	SUFFIX myexchangersoma
	USEION ca WRITE ica 
	USEION na WRITE ina
	RANGE ImaxNax, ica, ina , KnNacx, KcNacx, inacx
	GLOBAL dummy : prevent vectorization for use with CVODE
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
F = (faraday) (coulombs)
	R 	= (k-mole)	(joule/degC)
}

PARAMETER {
	ImaxNax = 0.511       (mA/cm2)     : <0,1e6>
	KnNacx   =  87.5     (mM)   <0,1e6>
	KcNacx   =  1.38     (mM)   <0,1e6>
}

ASSIGNED {
	celsius (degC)
	v (mV)
	ica (mA/cm2)
	ina (mA/cm2)
	cao (mM)
        cai (mM)
	nao (mM)
	nai (mM)
	inacx (mA/cm2) 
	dummy
}


BREAKPOINT { LOCAL q10, Kqa, KB, k

	k = R*(celsius + 273.14)/(F*1e-3)
		q10 = 3^((celsius - 37)/10 (degC))
	Kqa = exp(0.35*v/k)
	KB = exp( - 0.65*v/k)
			
	inacx = ImaxNax
	ina = -3*inacx   
	ica = 2*inacx          : 1*inacx
}
