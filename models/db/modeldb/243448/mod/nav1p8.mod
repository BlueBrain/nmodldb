TITLE TTX-R (Nav1.8) sodium current for bladder small DRG neuron soma model
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
        SUFFIX nav1p8
        USEION na READ ena WRITE ina
        RANGE gbar, g, ina
        RANGE minf, mtau, hinf, htau 
		THREADSAFE
}
 
PARAMETER {
        gbar = 0.0087177 (S/cm2)
}
 
STATE {
        m h
}
 
ASSIGNED {
		v		(mV)
		celsius (degC)
		ena		(mV)

		g		(S/cm2)
		ina		(mA/cm2)

		minf
		mtau	(ms)
		hinf
		htau	(ms)		
}
 

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gbar*m*m*m*h
	ina = g*(v - ena)
}
 

INITIAL {
	rates(v)
	
	m = minf
	h = hinf
}

DERIVATIVE states {  
        rates(v)
    
        m' = (minf-m)/mtau
		h' = (hinf-h)/htau 
}

PROCEDURE rates(v(mV)) {
                      
        LOCAL  alpha_m, beta_m, alpha_h, beta_h

UNITSOFF
		alpha_m = 7.21 -7.21/(1+exp((v-0.063)/7.86))
		beta_m = 7.4/(1+exp((v+53.06)/19.34))
		mtau = 1/(alpha_m+beta_m)			:Adapted from Han et al., 2015, Rat Nav .8 channel
		minf = 1/(1+exp((-11.4-v)/(8.5))) 	:Data Fit: Yoshimura et al., 1996

		alpha_h = 0.003 + 1.63/(1+exp((v+68.5)/10.01))
		beta_h = 0.81 -0.81/(1+exp((v-11.44)/13.12))
		htau = 1/(alpha_h+beta_h)			:Adapted from Han et al., 2015, Rat Nav 1.8 channel
		hinf = 1/(1+exp((v+24.2)/(5.6)))	:Data Fit: Yoshimura et al., 1996
}
UNITSON
