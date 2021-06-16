TITLE Cardiac fast sodium current
: Hodgkin - Huxley type sodium channel from Courtemanche et al Am J Physiol 1998 275:H301

NEURON {
	SUFFIX INaf
	USEION na READ ena WRITE ina
	RANGE gnabar, ina, Tauact, Tauinactf, Tauinacts
	GLOBAL minf, hinf, ninf, mtau, htau, ntau
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
        (mM) = (milli/liter)
        
}

PARAMETER {
	gnabar=.0156 (S/cm2) <0,1e9> 
	Tauact= 1 (ms) 
	Tauinactf=1 (ms) 
	Tauinacts=1 (ms)
   
               
}

STATE {
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
INITIAL {
	rates(v)
	m = minf
	h = hinf
        n = ninf   
	
}

BREAKPOINT {
	SOLVE states METHOD derivimplicit
	ina = gnabar*m*m*m*h*n*(v - ena)
}

DERIVATIVE states {
	rates(v)
	m' = (minf - m)/mtau
	h' = (hinf - h)/htau
        n' = (ninf - n)/ntau
}
UNITSOFF
FUNCTION alp(v(mV),i) (/ms) { LOCAL q10 :rest = -70  order m,h,j
	v = v 
	q10 = 3^((celsius - 37(degC))/10(degC))
	if (i==0) {
		alp = q10*0.32(/ms)*expMa(v *1(/mV) + 47.13, 47.13)
	}else if (i==1){
		alp = q10*0.135(/ms)*expMb(v *1(/mV) + 80, 6.8)
	}else if (i==2){
		alp = q10*1(/ms)*expMd(v *1(/mV), 79.23)
	}
}

FUNCTION bet(v(mV),i)(/ms) { LOCAL q10 :rest = -70  order m,h,n
	 v = v 
	q10 = 3^((celsius - 37(degC))/10(degC))
	if (i==0) {
		bet = q10*0.08(/ms)*exp(-v/11(mV))
	}else if (i==1){
		bet = q10*1(/ms)*expMc(v *1(/mV), 11.1)
	}else if (i==2){
		bet = q10*1(/ms)*expMe(v *1(/mV), 40.14)
	}
}

FUNCTION expMb(x,y) {
	if (v < -40) {
		expMb = (exp(-x/y))
	} else{
	expMb = 0
		
	}
}

FUNCTION expMc(x,y) {
	if (v < -40) {
		expMc = 3.56*exp(x*0.079) + 3.1e5*exp(0.35*x)
	}else{
		expMc = 1/(0.13*(1 + (exp(-(x + 10.66)/y))))
	}
}


FUNCTION expMd(x,y) {
	if (v < -40) {
		expMd = (-1.2740e5*exp(x*0.2444) - 3.474e-5*exp(-0.0439*x))*(x + 37.78)/(1 + exp(0.311*(x + y)))
	}else{
		expMd = 0
	}
}

FUNCTION expMe(x,y) {
	if (v < -40) {
		expMe = 0.1212*exp(-x*0.01052)/(1 + exp(-0.137*(x + y)))
	}else{
		expMe = 0.3*exp(-x*2.535e-7)/(1 + exp(-0.1*(x + 32)))
	}
}

FUNCTION expMa(x,y) {
	if (fabs(x/y) < 1e-6) {
		expMa = 10*(1 - x/y/2)
	}else{
		expMa = x/(1 - exp(-0.1*x))
	}
}


PROCEDURE rates(v) {LOCAL a, b
	:TABLE minf, hinf, ninf, mtau, htau, ntau DEPEND celsius FROM -100 TO 100 WITH 200
	a = alp(v,0)  b=bet(v,0)
	mtau = 1/(a + b)*Tauact
	minf = a/(a + b)
	a = alp(v,1)  b=bet(v,1)
	htau = 1/(a + b)*Tauinactf
	hinf = a/(a + b)
a = alp(v,2)  b=bet(v,2)
	ntau = 1/(a + b)*Tauinacts
	ninf = a/(a + b)
}
UNITSON
