TITLE K-D
: K-D current for prefrontal cortical neuron ------Yuguo Yu  2007

NEURON {
	SUFFIX kd
	USEION k READ ek WRITE ik
	RANGE gk, gbar, i
	GLOBAL minf, mtau, hinf, htau
}

PARAMETER {
	gbar = 0.1   	(S/cm2)	
								
	celsius
	ek = -100	(mV)          
	v 		(mV)
	vhalfm=-43  (mV)
	km=8
	vhalfh=-67  (mV) 
      kh=7.3
	q10=2.3
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

ASSIGNED {
	i  (mA/cm2)
	ik 		(mA/cm2)
	minf
	mtau (ms)	 	
	hinf
	htau (ms)	 	
}
 

STATE { m h}

BREAKPOINT {
        SOLVE states METHOD cnexp
       i = gbar*m*h*(v-ek)
	   ik=i
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
        minf=1-1/(1+exp((v-vhalfm)/km))
        hinf=1/(1+exp((v-vhalfh)/kh))

  	 mtau = 0.6
	 htau = 1500
}