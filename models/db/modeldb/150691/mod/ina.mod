TITLE Cardiac sodium current
: from BEELER & REUTER, J.Physiol, 1977

NEURON {
    THREADSAFE
	SUFFIX INa
	USEION na READ ena WRITE ina
	RANGE  gnabar, gnac, ina, Tauact, Tauinactf, Tauinacts
	GLOBAL minf, hinf, ninf, mtau, htau, ntau
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mM) = (milli/liter)
}

PARAMETER {
	gnabar=0.004 (S/cm2) <0,1e9> 
	gnac = 0.000003 (S/cm2)
	Tauact= 1 (ms) 
	Tauinactf=1 (ms) 
	Tauinacts=1 (ms)
}

STATE { : m h j
	m h n 
}

ASSIGNED {
	v (mV)
	celsius (degC) : 37
	ena (mV)
	ina (mA/cm2)
	minf hinf ninf
	mtau (ms)
	htau (ms)
	ntau (ms)  
}
LOCAL k
INITIAL { : m h j
	rates(v)
	m = minf
	h = hinf
	n = ninf
}

BREAKPOINT {
	SOLVE states METHOD derivimplicit
	ina = (gnabar*m*m*m*h*n + gnac)*(v - ena)
}

DERIVATIVE states {
	rates(v)
	m' = (minf - m)/mtau
	h' = (hinf - h)/htau
	n' = (ninf - n)/ntau
}
UNITSOFF
FUNCTION alp(v(mV),i) (/ms) { 
:LOCAL q10 :order m,h,j
	if (i==0) {
		alp = (- v - 47)/(exp(-0.1*(v + 47))-1)
	}else if (i==1){
		alp = 0.126*exp(-0.25*(v + 77))
	}else if (i==2){
		alp = 0.055*exp(-0.25*(v + 78))/(exp(-0.2*(v + 78))+1)
	}
}

FUNCTION bet(v(mV),i)(/ms) { 
:LOCAL q10 :order m,h,n
	if (i==0) {
		bet = 40*exp(-0.056*(v + 72))
	}else if (i==1){
		bet = 1.7/(exp(-0.082*(v + 22.5))+1)
	}else if (i==2){
		bet = 0.3/(exp(-0.1*(v + 32))+1)
	}
}

PROCEDURE rates(v) { : m h j
LOCAL a, b
	:TABLE minf, hinf, ninf, mtau, htau, ntau DEPEND celsius FROM -100 TO 100 WITH 200
	a = alp(v,0)  b=bet(v,0)
	mtau = 1/(a + b)
	minf = a/(a + b)
	a = alp(v,1)  b = bet(v,1)
	htau = 1/(a + b)
	hinf = a/(a + b)
	a = alp(v,2)  b=bet(v,2)
	ntau = 1/(a + b)
	ninf = a/(a + b)
}
UNITSON
