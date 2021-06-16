TITLE nap :Persistent sodium current

NEURON {
        SUFFIX nap
        USEION na READ ena WRITE ina
        RANGE  gna, ina, qna
}

UNITS {
	(mV) =	(millivolt)
	(mA) =	(milliamp)
        (molar) = (1/liter)
	(mM) =	(millimolar)
	PI   = (pi) (1)
	FARADAY	= 96485.309 (coul/mole)
	R = (k-mole) (joule/degC)
}

INDEPENDENT {t FROM 0 TO 1 WITH 100 (ms)}

PARAMETER {
	gna   = 0.0 (S/cm2)
	qfact = 3
}
 
ASSIGNED {
	ena	(mV)
        v 	(mV)
        ina	(mA/cm2)
        minf
	hinf	
	taum	(ms)
	tauh	(ms)
        diam    (um)
}

STATE { m h qna }
 
BREAKPOINT {
        SOLVE state METHOD cnexp
        ina = gna*m*h*(v-ena)
}
 
INITIAL {
	rates(v)
	m = minf
	h = hinf
	ina = gna*m*h*(v-ena)
	qna = 0
}

DERIVATIVE state { 
        rates(v)
        m' = (minf-m)/taum
        h' = (hinf-h)/(tauh/qfact)
	qna' = (-ina*diam*PI*(1e4)/FARADAY)/(diam*diam*PI/4)    
}
 
PROCEDURE rates(v (mV)) {  
	TABLE minf, hinf, taum, tauh
	FROM -200 TO 200 WITH 201

	minf = 1/(1+exp(-(v+48.7)/4.4))
	hinf = 1/(1+exp((v+48.8)/9.98))
	taum = 1/(0.091*(v+38)/(1-exp(-(v+38)/5))-0.062*(v+38)/(1-exp((v+38)/5)))
	
	if (v<=-60) {
		tauh = 3700+2000*1/(0.091*(v+22+38)/(1-exp(-(v+22+38)/5))-0.062*(v+22+38)/(1-exp((v+22+38)/5)))
	}
	else {
		tauh = 1200+8000*1/(0.091*(v+36+38)/(1-exp(-(v+36+38)/5))-0.062*(v+36+38)/(1-exp((v+36+38)/5)))
	}
}
