TITLE K-DR
: K-DR current for Retinal Ganglion Cell from Skaliora et al (1995)
: M.Migliore Dec. 2001

NEURON {
	SUFFIX kdrgc
	USEION k READ ek WRITE ik
	RANGE  gbar
	GLOBAL minf, mtau, hinf, htau
}

PARAMETER {
	gbar = 0.02   	(mho/cm2)	
								
	celsius
	ek		(mV)            : must be explicitly def. in hoc
	v 		(mV)
	a0m=0.04
	vhalfm=-30
	zetam=3
	gmm=0.8

	a0h=0.0001
	vhalfh=-58
	zetah=4
	gmh=0.6
	hmin=700
	
	q10=1
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
	hinf 		htau (ms)	 	
}
 

STATE { m h}

BREAKPOINT {
        SOLVE states METHOD cnexp
	ik = gbar*m*m*m*h* (v - ek)
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
        qt=q10^((celsius-24)/10)
        minf = 1/(1 + exp(-(v-0.2)/9.4))
	mtau = betm(v)/(a0m*(1+alpm(v)))

        hinf = 1/(1 + exp((v+52.3)/7.4))
	htau = beth(v)/(qt*a0h*(1+alph(v)))
	if (htau<hmin/qt) {htau=hmin/qt}
}

FUNCTION alpm(v(mV)) {
  alpm = exp(1.e-3*zetam*(v-vhalfm)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION betm(v(mV)) {
  betm = exp(1.e-3*zetam*gmm*(v-vhalfm)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION alph(v(mV)) {
  alph = exp(1.e-3*zetah*(v-vhalfh)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION beth(v(mV)) {
  beth = exp(1.e-3*zetah*gmh*(v-vhalfh)*9.648e4/(8.315*(273.16+celsius))) 
}
