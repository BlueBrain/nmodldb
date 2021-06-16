TITLE FH channel

COMMENT

    Frankenhaeuser - Huxley channels for Xenopus with Kirchoff's law for driving force

    Frankenhaeuser B, Huxley AF (1964) The action potential in the myelinated nerve fiber of Xenopus Laevis as computed on the basis of voltage clamp data. J Physiol 171:302-15

    Original: http://senselab.med.yale.edu/modeldb/ShowModel.asp?model=3507

    Modified by Christian Roessert: using Kirchoff's law for driving force

ENDCOMMENT

NEURON {
	SUFFIX fh
	USEION na READ ena WRITE ina
    NONSPECIFIC_CURRENT il
    RANGE gnabar, gl, el, gna, il, vsh
	GLOBAL minf, hinf, mtau, htau
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
}

PARAMETER {
	celsius (degC) : 14
	gnabar = 0 (S/cm2)	<0,1e9>
    gl = 0 (S/cm2)	   <0,1e9>
    el = 0 (mV)
    vsh = -70 (mV)

}
STATE {
	m h
}
ASSIGNED {
    v (mV)
    ena (mV)
    gna (S/cm2)
	ina (mA/cm2)
    il (mA/cm2)
    minf hinf 
	mtau (ms) htau (ms)
    inf[2]
	tau[2] (ms)
}

BREAKPOINT {
	LOCAL ghkna
	SOLVE states METHOD cnexp 
    
    gna = gnabar*m*m*h
	ina = gna*(v - ena)
    il = gl*(v - el)

}

INITIAL {
	mh(v*1(/mV))
	m = inf[0]
	h = inf[1]
}

? states
DERIVATIVE states {	: exact when v held constant
	mh(v*1(/mV))
	m' = (inf[0] - m)/tau[0]
	h' = (inf[1] - h)/tau[1]
}

UNITSOFF
FUNCTION alp(v(mV),i) { LOCAL a,b,c,q10 :rest = -70  order m,h
	v = v-vsh
	q10 = 3^((celsius - 20)/10)
	if (i==0) {
		a=.36 b=22. c=3.
		alp = q10*a*expM1(b - v, c)
	}else if (i==1){
		a=.1 b=-10. c=6.
		alp = q10*a*expM1(v - b, c)
	}
}

FUNCTION bet(v,i) { LOCAL a,b,c,q10 :rest = -70  order m,h
	v = v-vsh
	q10 = 3^((celsius - 20)/10)
	if (i==0) {
		a=.4  b= 13.  c=20.
		bet = q10*a*expM1(v - b, c)
	}else if (i==1){
		a=4.5  b= 45.  c=10.
		bet = q10*a/(exp((b - v)/c) + 1)
	}
}

FUNCTION expM1(x,y) {
	if (fabs(x/y) < 1e-6) {
		expM1 = y*(1 - x/y/2)
	}else{
		expM1 = x/(exp(x/y) - 1)
	}
}

PROCEDURE mh(v) {LOCAL a, b :rest = -70
	TABLE inf, tau DEPEND celsius FROM -100 TO 100 WITH 200
	FROM i=0 TO 1 {
		a = alp(v,i)  b=bet(v,i)
		tau[i] = 1/(a + b)
		inf[i] = a/(a + b)
	}
}
UNITSON
