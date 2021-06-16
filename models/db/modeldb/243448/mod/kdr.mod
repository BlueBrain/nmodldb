TITLE KDR Current for bladder small DRG neuron soma model
: Author: Darshan Mandge (darshanmandge@iitb.ac.in)
: Computational Neurophysiology Lab
: Indian Institute of Technology Bombay, India 

: For details refer: 
: A biophysically detailed computational model of bladder small DRG neuron soma 
: Darshan Mandge and Rohit Manchanda, PLOS Computational Biology (2018)


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
}

NEURON {
        SUFFIX kdr
        USEION k READ ek WRITE ik
        RANGE gbar, g, ik
        GLOBAL ninf, ntau
		THREADSAFE
}
 
PARAMETER {
         gbar = 0.002688 (S/cm2)
}
 
STATE {
        n
}
 
ASSIGNED {
        v (mV)
        ek (mV)

		g (S/cm2)
		ik (mA/cm2)

		ninf
		ntau (ms)
}
 

BREAKPOINT {
        SOLVE states METHOD cnexp
        g = gbar*n*n*n*n
		ik = g*(v - ek)     
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
        
        :"n" potassium activation system
        ninf = 1/(1+exp(-(v+35)/15.4)) :Sheets et al., 2007 supplement
		ntau = 7.14085 + (636.5557/(17.07686*sqrt(3.14/2)))*exp(-2*((v+19.99951)/17.07686)^2): by fitting data from Yoshimura et al., 2006
}
 
UNITSON
