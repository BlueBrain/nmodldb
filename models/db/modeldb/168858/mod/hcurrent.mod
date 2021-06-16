TITLE Anomalous rectifier current from Rhodes and Llinas, J Physiol 565:765-781, 2005.

COMMENT

        Implemented by Christina Weaver, 2007 (christina.weaver@mssm.edu),
        Including general kinetics parameters.

    from Rhodes and Llinas, J. Physiol. 565:765-781, (2005).

    minf = 1.0 / (1.0 + exp( (v-(-75))/5.5))
    tauM = 100 / ( exp(- (v-(-75))/11.0 ) + exp((v-(-75))/11) )
    if tauM < 5, tauM = 5.

    generalizing this:

    gH = gbarH * m * (v-erev)

    dm/dt = ( minf(V) - m ) / mtau(V)
    minf = 1.0 / (1.0 + exp( -2*a*(v-Vh))
    tauM = 1.0 / ( b * exp(a*(v-Vh) ) + b * exp(-a*(v-Vh)) )
    if tauM < tauMin, tauM = tauMin.

    where Vh = -75, a = -1 / 11, b = 1/100, tauMin = 5.

    This follows the general kinetics format also used by Av-Ron and Vidal, 1999.

    For comparison, Traub et al (2003) used the following equations:

    minf = 1 / (1+exp((v+75)/5.5))
    tauM = 1 / (exp(-14.6-0.086*v) + exp(-1.87+0.07*v))


ENDCOMMENT

UNITS { 
	(mV) = (millivolt) 
	(mA) = (milliamp) 
} 

NEURON { 
	SUFFIX ar
	NONSPECIFIC_CURRENT i
	RANGE gbar              
        RANGE Vh, tauMin
        RANGE a, b
}

PARAMETER { 
	gbar = 0.0 	(mho/cm2)
	erev = -55	(mV)  
        Vh = -75	(mV)
        a = -0.4090909	(/mV)
        b = 0.001	(1)
        tauMin = 5.0	(ms)
} 

ASSIGNED { 
	i 		(mA/cm2) 
	minf 		(1)
	mtau 		(ms) 
	v		(mV)
} 
STATE {
	m
}

INITIAL {
	rates(v)
	m = minf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	i = gbar * m * ( v - erev )
}

DERIVATIVE states {
	rates(v)
	m' = (minf-m)/mtau
}

UNITSOFF 

PROCEDURE rates(V (mV)) {
	minf = 1 / ( 1 + exp( -2 * a * ( V - Vh )) )
	mtau = 1 / b / (exp( a*(V-Vh)) + exp (-a*(V-Vh)) )

        if( mtau < tauMin ) { mtau = tauMin }
}

UNITSON

