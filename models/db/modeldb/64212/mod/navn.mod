
TITLE Na

NEURON {
	SUFFIX navn
	USEION na READ ena WRITE ina
	RANGE  gbar
	GLOBAL minf, mtau, hinf, htau
}

PARAMETER {
	gbar = 0.012  	(mho/cm2)	
								
	celsius
	ena		(mV)            : must be explicitly def. in hoc
	v 		(mV)
	
	a0m=0.3		
	vhalfm=-38.9	
	zetam=0.05	
	gmm=0.2		
	mmin=0.02	

	vm = -42	
	km = 8		

	a0h=0.03	
	vhalfh=-80	
	zetah=0.09	
	gmh=0.5		
	hmin=0.3	

	vh = -60	
	kh = 2		

	q10=2.3
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

ASSIGNED {
	ina 		(mA/cm2)
	minf 		mtau (ms)	 	
	hinf 		htau (ms)	 	
}
 

STATE { m h}

BREAKPOINT {
        SOLVE states METHOD cnexp
	ina = gbar*m^3*h* (v - ena)
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
	if (mtau<mmin/qt) {mtau=mmin/qt}

        hinf = 1/(1 + exp((v-vh)/kh))
	htau = beth(v)/(qt*a0h*(1+alph(v)))
	if (htau<hmin/qt) {htau=hmin/qt}
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
