TITLE Sodium Density Mechanism

UNITS {
	(mV) = (millivolts)
	(mA) = (milliamp)
	(S) = (siemens)
}

NEURON {
	SUFFIX na
	USEION na WRITE ina
	RANGE gmax, e
}

PARAMETER {
	gmax = 0.0009333 (S/cm2)
	e = 50 (millivolt)
}

ASSIGNED {
	v (mV)
	ena (mA)
	ina (mA/cm2)
	g (S/cm2)
}

STATE { m h }

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gmax * h * m^3
	ina = g * (v - e)
}

INITIAL {
	m = alpham(v)/(alpham(v) + betam(v))
    h = alphah(v)/(alphah(v) + betah(v))
}

DERIVATIVE states {
	m' = alpham(v) * (1-m) - betam(v) * m
	h' = alphah(v) * (1-h) - betah(v) * h
}

FUNCTION alpham(Vm (mV)) (/ms) {
	UNITSOFF
	alpham = 10 * exp(0.1 * (Vm + 34)) 
	UNITSON
}

FUNCTION betam(Vm (mV)) (/ms) {
	UNITSOFF
	betam = 10 * exp(-0.1 * (Vm + 34))
	UNITSON
}

FUNCTION alphah(Vm (mV)) (/ms) {
	UNITSOFF
	alphah = 0.05 * exp(-0.125 * (Vm + 29)) 	
	UNITSON
}

FUNCTION betah(Vm (mV)) (/ms) {
	UNITSOFF
	betah = 0.05 * exp(0.125 * (Vm + 29))
	UNITSON
}