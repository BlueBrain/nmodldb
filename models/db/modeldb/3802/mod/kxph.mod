TITLE kx
: K current for Salamander rod photoreceptor from Beech and Barnes (1989)
: M.Migliore Mar 2002

NEURON {
	SUFFIX kx
	USEION k READ ek WRITE ik
	RANGE  gbar
	GLOBAL minf, mtau
}

PARAMETER {
	gbar = 0.010   	(mho/cm2)	
	q10=3								
	ek		(mV)            : must be explicitly def. in hoc
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
	ik 		(mA/cm2)
	minf 		
	mtau (ms)	 	
}
 

STATE { m }

BREAKPOINT {
        SOLVE states METHOD cnexp
	ik = gbar*m * (v - ek)
} 

INITIAL {
	trates(v)
	m=minf  
}

DERIVATIVE states {   
        trates(v)      
        m' = (minf-m)/mtau
}

PROCEDURE trates(vm) {  
        LOCAL  a, b, qt
        qt=q10^((celsius-23)/10)
	a = trap0(vm,-45,0.23,7)		
	b = trap0(vm,-55,-1.2,-5.5)		
	minf = a/(a+b)
	mtau = 1.e3/((a+b)*qt)

}

FUNCTION trap0(v,th,a,q) {
	if (fabs(v-th) > 1e-6) {
	        trap0 = a * (v - th) / (1 - exp(-(v - th)/q))
	} else {
	        trap0 = a * q
 	}
}	

        

