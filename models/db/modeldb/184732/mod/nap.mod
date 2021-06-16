TITLE NAP Density Mechanism

UNITS {
	(mV) = (millivolts)
	(mA) = (milliamp)
	(S) = (siemens)
}

NEURON {
	SUFFIX nap
	USEION na WRITE ina
	RANGE gmax, e
}

PARAMETER {
	gmax = 0.000133333 (S/cm2)
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
	g = gmax * m * h
	ina = g * (v - ena)
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
	alpham = 0.25 * exp(0.1 * (Vm + 45)) 
	UNITSON
}

FUNCTION betam(Vm (mV)) (/ms) {
	UNITSOFF
	betam = 0.25 * exp(-0.083333 * (Vm + 45))
	UNITSON
}

FUNCTION alphah(Vm (mV)) (/ms) {
	UNITSOFF
	alphah = 5e-5 * exp(-0.083333 * (Vm + 53)) 	
	UNITSON
}

FUNCTION betah(Vm (mV)) (/ms) {
	UNITSOFF
	betah = 5e-5 * exp(0.083333 * (Vm + 53))
	UNITSON
}