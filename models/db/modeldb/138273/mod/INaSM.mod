TITLE Sodium channel
:  for spinal motoneuron

NEURON {
	SUFFIX INaSM
	USEION na READ ena WRITE ina
	RANGE gnamax, ina
}

UNITS {
	(mV)=(millivolt)
	(mA)=(milliamp)
}

PARAMETER {
	gnamax=0 (mho/cm2)
}

ASSIGNED {
	v (mV)
	ena (mV)
	ina (mA/cm2)   
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
    ina=gnamax*m*m*m*h*(v-ena)
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
UNITSOFF
    if (i==0){
          
            alfa=0.182*(v+40)/(1-exp(-(v+40)/9))
          
            }
    else if (i==1){
          
            alfa=0.024*(v+50)/(1-exp(-(v+50)/5))
          
            }
UNITSON
}

FUNCTION beta(v(mV),i) {
UNITSOFF
    if (i==0){
          
            beta=-0.124*(v+35)/(1-exp((v+35)/9))
          
            }
    else if (i==1){
          
            beta=-0.0091*(v+75)/(1-exp((v+75)/5))
          
            }
UNITSON
}

FUNCTION hinfin(v(mV)){
    UNITSOFF      
            hinfin=1/(1+exp((v+65)/6.2))
UNITSON
          
}

PROCEDURE settables(v (mV)) {LOCAL a,b,c
    TABLE taum, minf, tauh, hinf FROM -100 TO 100 WITH 2000
UNITSOFF
    a=alfa(v,0) b=beta(v,0)
    taum=1/(a+b)
    minf=a/(a+b)
    a=alfa(v,1) b=beta(v,1) c=hinfin(v)
    tauh=1/(a+b)
    hinf=c
UNITSON
}
UNITSON