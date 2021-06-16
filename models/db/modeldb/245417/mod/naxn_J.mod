TITLE nax Engel and Jonas 2005
: Na current for Mossy Fiber Boutons.
: Na current (M.Migliore Jul. 1997 and Apr.2002) modified by E.Giacalone following Engel and Jonas 2005 sodium channel kinetics

NEURON {
	SUFFIX naxj
	USEION na READ ena WRITE ina
	RANGE  gbar, sh
	GLOBAL minf, hinf, mtau, htau :thinf, qinf
}

PARAMETER {
	sh   = 0	(mV)
	gbar = 0.010   	(mho/cm2)	
								
	tha  = -105.023 (mV)		: v 1/2 for act	
	qa   = 17.7094	(mV)		: act slope (4.5)		
	Ra   = 93.8285	(/ms)		: open (v)		
	Rb   = 0.168396	(/ms)		: close (v)		
	qb   = 23.2707	(mV)		: act slope (4.5)

		
	thi  = 17.6769	(mV)		: v 1/2 for inact 	
	qd   = 18.706	(mV)	        : inact tau slope
	qg   = 13.3097  (mV)
	mmin=0.002	
	hmin=0.05			
	q10=2
	Rg   = 6.62694 	(/ms)		: inact recov (v) 	
	Rd   = 0.000354	(/ms)		: inact (v)	

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
	thegna		(mho/cm2)
	minf 		hinf 		
	mtau (ms)	htau (ms) 	
}
 
STATE { m h}

BREAKPOINT {
        SOLVE states METHOD cnexp
        thegna = gbar*m*m*m*h
	ina = thegna * (v - ena)
} 

INITIAL {
	trates(v,sh)
	m=minf  
	h=hinf
}

DERIVATIVE states {   
        trates(v,sh)      
        m' = (minf-m)/mtau
        h' = (hinf-h)/htau
}


PROCEDURE trates(vm,sh2) {  
        LOCAL  a, b, qt
        qt=q10^((celsius-24)/10)
	a = trap0(vm,-tha-sh2,Ra,qa) 
	b = Rb * exp(-vm/qb)
	mtau = 1/(a+b)/qt
        if (mtau<mmin) {mtau=mmin}
	minf = a/(a+b)

	a = Rd * exp(-vm/qd)
	b = Rg/ (exp(-(v+thi)/qg)+1)
	htau =  1/(a+b)/qt
        if (htau<hmin) {htau=hmin}
	hinf = a/(a+b)
}

FUNCTION trap0(v,th,a,q) {
	if (fabs(v-th) > 1e-6) {
	        trap0 = a * (v - th) / (1 - exp(-(v - th)/q))
	} else {
	        trap0 = a * q
 	}
}	

