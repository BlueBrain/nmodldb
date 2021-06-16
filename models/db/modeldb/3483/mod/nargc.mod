TITLE na3
: Na current for Retinal Ganglion Cell from Benison et al (2001)
: M.Migliore Nov. 2001

NEURON {
	SUFFIX nargc
	USEION na READ ena WRITE ina
	RANGE  gbar
	GLOBAL minf, hinf, mtau, htau
}

PARAMETER {
	gbar = 0.010   	(mho/cm2)	
								
	ena		(mV)            : must be explicitly def. in hoc
	celsius
	v 		(mV)
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

ASSIGNED {
	ina 		(mA/cm2)
	minf 		hinf 		
	mtau (ms)	htau (ms) 	
}
 

STATE { m h}

BREAKPOINT {
        SOLVE states METHOD cnexp
	ina = gbar*m*m*m*h * (v - ena)
} 

INITIAL {
	trates(v)
	m=minf  
	h=hinf
}

DERIVATIVE states {   
        trates(v)      
        m' = (minf-m)/mtau
        h' = (hinf-h)/htau
}

PROCEDURE trates(vm) {  
        LOCAL  a, b

	a = trap0(vm,-29,0.5,0.18)
	b = 6*exp(-(vm+45)/15)
	minf = a/(a+b)
	mtau = 1/(a+b)

	a = 0.15*exp(-(vm+47)/20)
	b = 2.8/(1+exp(-0.1*(vm+20)))
	hinf = a/(a+b)
	htau =  1/(a+b)
}

FUNCTION trap0(v,th,a,q) {
	if (fabs(v-th) > 1e-6) {
	        trap0 = a * (v - th) / (1 - exp(-(v - th)*q))
	} else {
	        trap0 = a / q
 	}
}	

        

