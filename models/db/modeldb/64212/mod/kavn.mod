TITLE K-A
: K-A current for Mitral Cells from Wang et al (1996)
: M.Migliore Jan. 2002
: adapted for vn neurons M.Migliore 2005

NEURON {
	SUFFIX kavn
	USEION k READ ek WRITE ik
	RANGE  gbar
	GLOBAL minf, mtau, hinf, htau
}

PARAMETER {
	gbar = 0.00215   	(mho/cm2)	
								
	celsius
	ek		(mV)            : must be explicitly def. in hoc
	v 		(mV)

	a0m=0.0035	
	vhalfm=-75	
	zetam=0.2	
	gmm=0.82 	

	vm=-21		
	km=12.8		

	a0h=0.002	
	vhalfh=-70	
	zetah=0.05	
	gmh=0.95	

	vh=0		
	kh=60		

	q10=3
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
	ik = gbar*m^1.5*h*(v - ek)
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
        qt=q10^((celsius-22)/10)
        minf = 1/(1 + exp(-(v-vm)/km))
	mtau = betm(v)/(qt*a0m*(1+alpm(v)))

        hinf = 1/(1 + exp((v-vh)/kh))
	htau = beth(v)/(qt*a0h*(1+alph(v)))
}

FUNCTION alpm(v(mV)) {
  alpm = exp(zetam*(v-vhalfm)) 
}

FUNCTION betm(v(mV)) {
  betm = exp(zetam*gmm*(v-vhalfm)) 
}

FUNCTION alph(v(mV)) {
  alph = exp(zetah*(v-vhalfh)) 
}

FUNCTION beth(v(mV)) {
  beth = exp(zetah*gmh*(v-vhalfh)) 
}
