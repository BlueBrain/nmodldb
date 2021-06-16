TITLE HH sodium channel
: Hodgkin - Huxley squid sodium channel
: file update to provide temperature dependence 1/17/2006

NEURON {
	SUFFIX na_sej
	USEION na READ ena WRITE ina
	RANGE gnasejbar, ina, gna

}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
        v (mV)
        dt (ms)
	gnasejbar=.086 (mho/cm2) <0,1e9>
        ena = 55 (mV)
}

STATE {
	m h
}

ASSIGNED {
	ina (mA/cm2)
	minf hinf
	mtau (ms)
	htau (ms)
        gna (mho/cm2)
}

INITIAL {
	rate(v)
	m = minf
	h = hinf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gna = gnasejbar*m*m*m*h
        ina = gna*v*(15/135 - exp(-v/25.4))/(1-exp(-v/25.4))
}

DERIVATIVE states {
	rate(v)
	m' = (minf - m)/mtau
	h' = (hinf - h)/htau
}

UNITSOFF

FUNCTION malf(v(mV))(/ms){ LOCAL va
	va = v + 35  
	if (fabs(va)<1e-04) {
		malf = -0.182*(-9 + 0.5*va)
	}else{
		malf = 0.182*(v+35)/(1-exp(-(v+35)/9)) 
	}
}

FUNCTION mbet(v(mV))(/ms) { LOCAL vb
	vb = v + 35
	if (fabs(vb)<1e-04) {
		mbet = 0.124*(9+vb*0.5)
	}else{
		mbet = 0.124*(v+35)/(-1+exp((v+35)/9))
	}
}

FUNCTION half(v(mV))(/ms) { LOCAL vc
	vc = v + 48 
	if (fabs(vc)<1e-04) {
		half = -0.025*(4 + 0.5*vc)
	}else{
		half = 0.025*(v+48)/(1-exp(-(v+48)/4)) 
	}
}

FUNCTION hbet(v(mV))(/ms) { LOCAL vd
	vd = v + 72
	if (fabs(vd)<1e-04) {
		hbet = 0.0091*(5+vd*0.5)
	}else{
		hbet = 0.0091*(v+72)/(-1+exp((v+72)/5))
	}
}



PROCEDURE rate(v(mV)) {LOCAL q10, q11, msum, hsum, ma, mb, ha, hb
	TABLE minf, hinf, mtau, htau DEPEND celsius FROM -100 TO 100 WITH 200

	q10 = (2.8)^((celsius - 23)/10)
        q11 = (2.4)^((celsius - 23)/10) 

	ma=malf(v) mb=mbet(v) ha=half(v) hb=hbet(v) 
	msum = ma + mb
        minf = ma/msum
        mtau = 1/(q10*msum)
        
        hsum = ha + hb
        hinf = 1/(1+exp((v+65)/6.2))
        htau = 1/(q11*hsum)


}

UNITSON
