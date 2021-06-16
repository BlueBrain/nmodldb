TITLE KA
: K-A current for hippocampal interneurons from Lien et al (2002)
: M.Migliore Jan. 2003

NEURON {
	SUFFIX ka
	USEION k READ ek WRITE ik
	RANGE  gbar
	GLOBAL minf, hinf, htau, mtau
}

PARAMETER {
	gbar = 0.0002   	(mho/cm2)	
								
	celsius
	ek		(mV)            : must be explicitly def. in hoc
	v 		(mV)
	a0h=0.17
	vhalfh=-105
	q10=3
	hmin=5
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

ASSIGNED {
	ik 		(mA/cm2)
	minf 		mtau (ms)
	hinf	 	htau (ms)
}
 

STATE { m h}

BREAKPOINT {
        SOLVE states METHOD cnexp
	ik = gbar*m*h*(v - ek)
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

PROCEDURE trates(v) {  
	LOCAL qt
        qt=q10^((celsius-23)/10)
        minf = (1/(1 + exp(-(v+41.4)/26.6)))^4
	mtau=0.5/qt
        hinf = 1/(1 + exp((v+78.5)/6))
	htau = a0h*(v-vhalfh)/qt
	if (htau<hmin/qt) {htau=hmin/qt}
}

