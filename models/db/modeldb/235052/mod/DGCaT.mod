TITLE T-type calcium channel based originally on deep cerebellar nucleus (DCN) neuron


NEURON { 
	SUFFIX DGCaT 
	USEION ca WRITE ica
	RANGE m, h, cali, taum, tauh, mshift, hshift, ic, base, gcat
	GLOBAL qdeltat
} 
 
UNITS { 
	(mA) = (milliamp) 
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
} 
 
PARAMETER { 
    qdeltat = 1
    gcat = 1 (S/cm2)
    mshift (mV)
    hshift (mV)

    sep (mV)
  
    cali = 0.00007 (mM)
    calo = 2  (mM)   
    
    base = 9.3 ()  
    
    ic (mA/cm2)
} 

ASSIGNED {
    v (mV)



	ica (mA/cm2) 
	minf
	hinf
	taum (ms) 
	tauh (ms) 
	celsius (degC)
	T (kelvin)
    A (1)
} 
 
STATE {
	m
    h
} 

INITIAL { 
    T = 273.15 + celsius
    rate(v)
    m = minf 
	h = hinf
} 
 
BREAKPOINT { 
    SOLVE states METHOD cnexp 
  
    ica = gcat * m*m * h * (v-115)
    ic = ica
} 
 
DERIVATIVE states { 
	rate(v) 
	m' = (minf - m)/taum 
	h' = (hinf - h)/tauh 
} 

PROCEDURE rate(v(mV)) {
	TABLE minf, taum, hinf, tauh  FROM -150 TO 100 WITH 300 
	minf = 1 / (1 + exp(-.16*(v + 55 - mshift)))
	taum = 0 + (0.2*exp(-.049*(v-mshift)))
	
	hinf = 1 / (1 + exp(0.25*(v + 80+hshift)))    : 0.25 scalar
   	tauh = base + (0.045*exp(-0.095*(v+hshift)))
}

