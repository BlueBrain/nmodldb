TITLE ECellChann  
: Cell model from  Olufsen et al. (2003)
: Programmed by Adriano Tort, CBD, BU, 2008

UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
	SUFFIX Ecellchann
	NONSPECIFIC_CURRENT ina
	NONSPECIFIC_CURRENT ik
:	NONSPECIFIC_CURRENT celldrive,celldrivenoise
	NONSPECIFIC_CURRENT il
:	GLOBAL gna, gk, ena, ek, gl, el
	RANGE gna, gk, ena, ek, gl, el
	RANGE drive,drivenoise 
}
 
PARAMETER {
        v				(mV)
        celsius			(degC)
        gna= 0.1		(mho/cm2)
        ena= 50			(mV)
        gk= 0.080		(mho/cm2)
        ek= -100		(mV)
		gl= 0.0001		(mho/cm2)
        el= -67			(mV)
:		drive=0.0004	(mA/cm2)
:		drivenoise =0	(mA/cm2)
}
 
STATE {
	m
	n
	h 
}
 
ASSIGNED {
	celldrive		(mA/cm2)
	celldrivenoise	(mA/cm2)
	
	ina		(mA/cm2) 
 	minf
	mtau    (ms)
	hinf
	htau	(ms)
	
	ik		(mA/cm2)
	ninf
	ntau	(ms)
	il		(mA/cm2) 
}
 
BREAKPOINT {
        SOLVE states METHOD cnexp
:		celldrive = -drive
:		celldrivenoise = -drivenoise
        
		ik=gk*n*n*n*n*(v-ek)
		ina=gna*minf*minf*minf*h*(v-ena)
		il=gl*(v-el)      
}
 
DERIVATIVE states { 
       rates(v)
:	   n'= (ninf- n)/ ntau 
:	   h'= (hinf- h)/ htau 
	   n'= (ninf- n)/ ntau /4
	   h'= (hinf- h)/ htau /4
}


INITIAL { 
	rates(v)
	n=ninf
	h=hinf
}

PROCEDURE rates(v (mV)) {
LOCAL alpha, beta
UNITSOFF 

	alpha = 0.32*(v+54)/(1-exp(-(v+54)/4))
	beta = 0.28*(v+27)/(exp((v+27)/5)-1)
	minf=alpha/(alpha+beta)
	
	
	alpha = 0.128*exp(-(v+50)/18)
  	beta =  4/(1 + exp(-(v+27)/5))
	htau = 1/(alpha+beta)
	hinf = alpha*htau	
	
	
	alpha = 0.032*(v+52)/(1-exp(-(v+52)/5))
	beta = 0.5*exp(-(v+57)/40)
	ntau=1/(alpha+beta)
	ninf=alpha*ntau	
	
UNITSON
}
