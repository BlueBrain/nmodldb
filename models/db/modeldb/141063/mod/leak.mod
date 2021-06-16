TITLE passive membrane channel
: Added temperature-dependency
: Written by Patricio Orio

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
	(S) = (siemens)
}

NEURON {
	SUFFIX leak
	NONSPECIFIC_CURRENT i
	RANGE g, e
}

PARAMETER {
	g = .0001	(S/cm2)	<0,1e9>
	e = -70	(mV)
}

ASSIGNED {
	v (mV)
	i (mA/cm2)
	celsius (degC)
	rho
}

BREAKPOINT {
	rho = 1.3 ^ ((celsius - 25 (degC))/10 (degC))
	i = rho * g * (v - e)
}
