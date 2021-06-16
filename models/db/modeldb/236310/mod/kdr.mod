TITLE KDR (4ap resistant, persistent) current for Evans 2012

COMMENT


Recorded at 22C - corrected to 35C with qfact 3

Jason Moyer 2004 - jtmoyer@seas.upenn.edu

ENDCOMMENT

UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (S)  = (siemens)
}
 
NEURON {
        SUFFIX kdr
        USEION k READ ek WRITE ik
        RANGE  gkbar, ik
}
 
PARAMETER {
	gkbar   =   0.0191 (S/cm2)

	avhalf = -13		(mV)	: Nisenbaum 1996, Fig 6C
	aslope = -9.09		(mV)	: Nisenbaum 1996, Fig 6C
	ashift = 13.605			(mV)
	bvhalf = -13		(mV)
	bslope = -12.5		(mV)
	bshift = 0.33739			(mV)


 	a = 0.7				: matched to Nisenbaum 1996, figure 9A (with qfact = 1)
 	qfact = 0.36 
}
 
STATE { m }
 
ASSIGNED {
	ek				(mV)
        v 				(mV)
        ik 				(mA/cm2)
        gk				(S/cm2)
        minf 
        alpham
        betam
        taum		(ms)
    }
 
BREAKPOINT {
        SOLVE state METHOD cnexp
        gk = gkbar * m  
        ik = gk * ( v - ek )
}
 

 
INITIAL {
	rates(v)
	
	m = minf
}


DERIVATIVE state { 
        rates(v)
        m' = (minf - m) / (taum/qfact)

}
 
PROCEDURE rates(v (mV)) {  
			alpham=exp((v-avhalf-ashift)/aslope)
			betam=exp((v-bvhalf-bshift)/bslope)
			taum=(0.01*50*betam/(1+alpham))/qfact
			minf=1/(1+alpham)		
}
 
 
