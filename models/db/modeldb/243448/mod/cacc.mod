TITLE CACC current for bladder small DRG neuron soma model

: Author: Darshan Mandge (darshanmandge@iitb.ac.in)
: Computational Neurophysiology Lab
: Indian Institute of Technology Bombay, India 

: For details refer: 
: A biophysically detailed computational model of bladder small DRG neuron soma 
: Darshan Mandge and Rohit Manchanda, PLOS Computational Biology (2018)

NEURON {
	SUFFIX cacc
	USEION cl READ ecl WRITE icl VALENCE -1
	USEION caip3r READ caip3ri VALENCE 2
	RANGE g, icl  
	RANGE gbar, ninf, ntau, vhalf, sf1, EC50, hc
}

UNITS {
	(molar) = (1/liter)
	(mM)	= (millimolar)
	(uM)	= (micromolar)
	(S)  	= (siemens)
	(mA) 	= (milliamp)
	(mV) 	= (millivolt)
}

PARAMETER {
        gbar = 1e-6 (S/cm2)
}

ASSIGNED {
        v       (mV)
        icl		(mA/cm2)
		g		(mho/cm2)
		ninf
		ntau 	(ms)
		caip3ri	(mM)
		celsius	(degC)
		ecl 	(mV)
		EC50 	(mM)
		hc 		(1)
}

STATE {
        n
}
 
BREAKPOINT {
        SOLVE states METHOD cnexp
		g = gbar*n
		icl = g*(v - ecl)
}
 
INITIAL {
		rates(v)
		n = ninf
}

DERIVATIVE states {  
        rates(v)
        n' = (ninf-n)/ntau
}
 
PROCEDURE rates(v(mV)) {
		 
		 UNITSOFF
        :"n" CACC activation 		
		hc = -0.3126*exp(-v/81.02) + 2.086
		EC50 = 0.39175*exp(-v/38.307) + 0.468
		
		if(v<-100){ 
			hc = 1.012
			EC50 = 5.798
		}
		
		ninf = (1/(1+(EC50*1e-3/caip3ri)^hc))
		ntau = 1
}
 
UNITSON 