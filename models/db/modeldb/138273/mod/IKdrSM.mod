TITLE Potassium delayed rectifier channel
:  for spinal motoneuron

NEURON {
	SUFFIX IKdrSM
	USEION k READ ek WRITE ik
	RANGE gkdrmax, ik
}

UNITS {
	(mV)=(millivolt)
	(mA)=(milliamp)
}

PARAMETER {
	gkdrmax=0 (mho/cm2)
}

ASSIGNED {
	v (mV)
	ek (mV)
	ik (mA/cm2)   
	minf
	taum (ms)	
}

STATE {
	m
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ik=gkdrmax*m*m*m*m*(v-ek)
}

INITIAL {
    settables(v)
	m=minf
}

DERIVATIVE states {
    settables(v)
    m'=(minf-m)/taum
}

FUNCTION alfa(v(mV)) {
 UNITSOFF           
            alfa=0.0075*(v+30)/(1-exp(-(v+30)/10))
UNITSON            
}

FUNCTION beta(v(mV)) {
UNITSOFF            
            beta=0.1*exp(-(v+46)/31)
UNITSON            
}



PROCEDURE settables(v (mV)) {LOCAL a,b,c
  UNITSOFF
    TABLE taum, minf FROM -100 TO 100 WITH 2000
    a=alfa(v) b=beta(v)
    taum=1/(a+b)
    minf=a/(a+b)
  UNITSON
}
UNITSON