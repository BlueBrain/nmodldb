TITLE Potassium fast channel
:  for spinal motoneuron

NEURON {
	SUFFIX IKaSM
	USEION k READ ek WRITE ik
	RANGE gkamax, ik
}

UNITS {
	(mV)=(millivolt)
	(mA)=(milliamp)
}

PARAMETER {
	gkamax=0 (mho/cm2) 
}

ASSIGNED {
	v (mV)
	ek (mV)
	ik (mA/cm2)   
	minf
	taum (ms)
	hinf
	tauh (ms)
}

STATE {
	m h
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ik=gkamax*m*m*m*m*h*(v-ek)
}

INITIAL {
    settables(v)
	m=minf
	h=hinf
}

DERIVATIVE states {
    settables(v)
    m'=(minf-m)/taum
    h'=(hinf-h)/tauh
}

FUNCTION alfa(v(mV),i) {
    if (i==0){
            UNITSOFF
            alfa=0.032*(v+64)/(1-exp(-(v+64)/6))
            UNITSON
            }
    else if (i==1){
            UNITSOFF
            alfa=0.05/(1+exp((v+86)/10))
            UNITSON
            }
}

FUNCTION beta(v(mV),i) {
    if (i==0){
            UNITSOFF
            beta=0.203*exp(-(v+40)/24)
            UNITSON
            }
    else if (i==1){
            UNITSOFF
            beta=0.05/(1+exp(-(v+86)/10))
            UNITSON
            }
}


PROCEDURE settables(v (mV)) {LOCAL a,b
    UNITSOFF
    TABLE taum, minf, tauh, hinf FROM -100 TO 100 WITH 2000
    a=alfa(v,0) b=beta(v,0)
    taum=1/(a+b)
    minf=a/(a+b)
    a=alfa(v,1) b=beta(v,1)
    tauh=1/(a+b)
    hinf=a/(a+b)
    UNITSON
}
UNITSON