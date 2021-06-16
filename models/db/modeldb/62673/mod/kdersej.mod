TITLE HH sodium channel
: Hodgkin - Huxley squid sodium channel
: file updated to provide temperature dependence 1/17/2006

NEURON {
	SUFFIX kder_sej
	USEION k READ ek WRITE ik
	RANGE gkdersejbar, ik, gkder

}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
        v (mV)
        dt (ms)
	gkdersejbar=.086 (mho/cm2) <0,1e9>
        ek = -104 (mV)
}

STATE {
	m
}

ASSIGNED {
	ik (mA/cm2)
	minf hinf
	mtau (ms)
        gkder (mho/cm2)
}

INITIAL {
	rate(v)
	m = minf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gkder = gkdersejbar*m*m*m*m
        ik = gkder*(v - ek)
}

DERIVATIVE states {
	rate(v)
	m' = (minf - m)/mtau

}

UNITSOFF

FUNCTION malf(v(mV))(/ms){ LOCAL va
	va = v + 20  
	if (fabs(va)<1e-04) {
		malf = -0.02*(-9 + 0.5*va)
	}else{
		malf = 0.02*(v+20)/(1-exp(-(v+20)/9)) 
	}
}

FUNCTION mbet(v(mV))(/ms) { LOCAL vb
	vb = v + 20
	if (fabs(vb)<1e-04) {
		mbet = 0.002*(9+vb*0.5)
	}else{
		mbet = 0.002*(v+20)/(-1+exp((v+20)/9))
	}
}




PROCEDURE rate(v(mV)) {LOCAL q10, msum, ma, mb
	TABLE minf, mtau DEPEND celsius FROM -100 TO 100 WITH 200

        q10 = (2.8)^((celsius - 23)/10)
	ma=malf(v+1) mb=mbet(v+1) 
	msum = ma + mb
        minf = ma/msum
        mtau = 1/(q10*msum)



}

UNITSON
