TITLE Slowly recovering Na channel  
: Na channel in Large Cell of Corpus Glomerulosum
: Reconstruction of physiological data    :Kinetic parameters updated 1 May 2001 
 
NEURON {
	SUFFIX Na_cglc
	USEION na READ ena WRITE ina
	RANGE gnabar, ina
	GLOBAL minf, hinf, mtau, htau
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gnabar=.120 (mho/cm2) <0,1e9>
	
}

STATE {
	m h
}

ASSIGNED {
	v (mV)
	celsius (degC) : 6.3
	ena (mV)
	ina (mA/cm2)
	minf hinf
	mtau (ms)
	htau (ms)
}

INITIAL {
	rates(v)
	m = minf
	h = hinf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ina = gnabar*m*m*m*h*(v - ena)
}

DERIVATIVE states {
	rates(v)
	m' = (minf - m)/mtau
	h' = (hinf - h)/htau
}

FUNCTION alp(v(mV),i) (/ms) { LOCAL a,b,c,q10 :rest = -70  order m,h
	
	q10 = 1		 :temp dependence not incorporated
	if (i==0) {	 :m gate	
		alp = q10*(.03505(/ms)*1(/mV)*(v+42.29529)+1(/ms)*(0.0012285*(1(/mV)*v+42.29529)^2+0.005)^.5)
	}else if (i==1){ :h gate	
		alp = q10*.000187291(/ms)*exp(-v/20.8036(mV))
	}
}

FUNCTION bet(v(mV),i)(/ms) { LOCAL a,b,c,q10 :rest = -70  order m,h
	
	q10 = 1
	if (i==0) {
		bet = q10*.403703(/ms)*(1 - 1/(1 +exp((-1(/mV)*v - 44.6804)/10.0433)))
	}else if (i==1){
		bet = q10*.424283(/ms)/(1 + exp((-1(/mV)*v - 38.85) / 5.74763 ))
	}
}

FUNCTION expM1(x,y) {
	if (fabs(x/y) < 1e-6) {
		expM1 = y*(1 - x/y/2)
	}else{
		expM1 = x/(exp(x/y) - 1)
	}
}

PROCEDURE rates(v(mV)) {LOCAL a, b
	TABLE minf, hinf, mtau, htau DEPEND celsius FROM -100 TO 100 WITH 200
	a = alp(v,0)  b=bet(v,0)
	mtau = 1/(a + b)
	minf = a/(a + b)
	a = alp(v,1)  b=bet(v,1)
	htau = 1/(a + b)
	hinf = a/(a + b)
}
