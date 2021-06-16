TITLE Icellchann
: Cell model from Wang and Buzsaki, J Neurosci 1996
: Programmed by Adriano Tort, CBD, BU, 2008
: With modification by Ruben Tikidji-Hamburyan rth@nisms.krinc.ru
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
	SUFFIX BSKCch
	USEION na READ ena WRITE ina
	USEION k  READ ek  WRITE ik
	RANGE  gnabar, gkbar
	: , v0

}
 
PARAMETER {
	v				(mV)
	celsius			(degC)
	gna= 0.035		(mho/cm2)
	ena				(mV)
	gk= 0.009		(mho/cm2)
	ek				(mV)
:		gl= 0.0001		(mho/cm2)
:       el= -65			(mV)
:	v0=-65			(mV)
		}
 
STATE {
	m
	n
	h 
}
 
ASSIGNED {
	ina		(mA/cm2) 
 	minf
	mtau   	(ms)
	hinf
	htau    (ms)
	
	ik		(mA/cm2)
	ninf
	ntau	(ms)	
}
 
BREAKPOINT {
        SOLVE states METHOD cnexp
        ina=gna*minf*minf*minf*h*(v-ena)
		ik=gk*n*n*n*n*(v-ek)
}
 
DERIVATIVE states { 
       rates(v)
	   h'= 5*(hinf- h)/ htau
	   n'= 5*(ninf- n)/ ntau 
}


INITIAL { 
	rates(v)
	:DB>>
	:printf("BSKCch init v=%g\n",v)
	:<<DB
	n=0
	h=1
	}

PROCEDURE rates(v (mV)) {
LOCAL alpha, beta
UNITSOFF 

	alpha = -0.1*(v+35)/(exp(-(v+35)/10)-1)
	beta = 4*exp(-(v+60)/18)
	minf=alpha/(alpha+beta)
	
	alpha = 0.07*exp(-(v+58)/20)
	beta = 1/(exp(-0.1*(v+28))+1)
	htau=1/(alpha+beta)
	hinf=alpha*htau
	
	alpha = 0.01*(v+34)/(1-exp(-0.1*(v+34)))
	beta = 0.125*exp(-(v+44)/80)
	ntau=1/(alpha+beta)
	ninf=alpha*ntau	
	
UNITSON
}
