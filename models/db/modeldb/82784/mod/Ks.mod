: Slowly inactivating K+ channel
: from Durstewitz & Gabriel (2006), Cerebral Cortex

NEURON {
	SUFFIX Ks
	USEION k READ ki, ko WRITE ik
	RANGE gKsbar, gk
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
        (mM) = (milli/liter)
}

PARAMETER {
	gKsbar= 0.00014 (mho/cm2) <0,1e9>
}

STATE {
	a 
        b 
}

ASSIGNED {
        v  (mV)
	ik (mA/cm2)
	ainf
	binf
	atau (ms)
	btau (ms)
	gk (mho/cm2)
	ek  (mV)
	ki (mM)
	ko  (mM)
}

INITIAL {
	rate(v)
	a = ainf
	b = binf
}

BREAKPOINT {
	SOLVE states METHOD derivimplicit
	gk = gKsbar*a*b
	ek = 25*log(ko/ki)
	ik = gk*(v-ek)
}

DERIVATIVE states {
	rate(v)
	a' = (ainf-a)/atau
	b' = (binf-b)/btau
}

UNITSOFF

PROCEDURE rate(v(mV)) {
	ainf = 1/(1+exp(-(v+34)/6.5))
	atau = 6
	binf = 1/(1+exp((v+65)/6.6))
	btau = 200+220/(1+exp(-(v+71.6)/6.85))
}
	
UNITSON
