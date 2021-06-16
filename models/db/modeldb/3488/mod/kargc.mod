TITLE K-A
: K-A current for Retinal Ganglion Cell from Benison et al (2001)
: M.Migliore Nov. 2001

NEURON {
	SUFFIX kargc
	USEION k READ ek WRITE ik
	RANGE  gbar
	GLOBAL minf, hinf, mtau, htau
}

PARAMETER {
	gbar = 0.010   	(mho/cm2)	
								
	ek		(mV)            : must be explicitly def. in hoc
	v 		(mV)
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

ASSIGNED {
	ik 		(mA/cm2)
	minf 		hinf 		
	mtau (ms)	htau (ms) 	
}
 

STATE { m h}

BREAKPOINT {
        SOLVE states METHOD cnexp
	ik = gbar*m*m*m*h * (v - ek)
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

	a = trap0(vm,-15,0.02,0.12)
	b = 0.05*exp(-(vm+1)/30)
	minf = a/(a+b)
	mtau = 1/(a+b)

	hinf = 1/(1+exp((vm+62)/6.35))
	htau =  25
}

FUNCTION trap0(v,th,a,q) {
	if (fabs(v-th) > 1e-6) {
	        trap0 = a * (v - th) / (1 - exp(-(v - th)*q))
	} else {
	        trap0 = a / q
 	}
}	

        

