TITLE Sodium transient current for Yu et al. 2008

COMMENT

	Implemented by Erin Munro
ENDCOMMENT

INDEPENDENT { t FROM 0 TO 1 WITH 1 (ms) }

UNITS { 
	(mV) = (millivolt) 
	(mA) = (milliamp) 
} 
NEURON { 
	SUFFIX nafYu
	USEION na READ ena WRITE ina
	RANGE gbar, ina, m, h, df, am, bm, ah, bh, hinf
}
PARAMETER { 
	gbar = 0.0 	   (mho/cm2)
	v (mV) ena 		   (mV)  
} 
ASSIGNED { 
	ina 		   (mA/cm2) 
	am (1/ms) bm 	   (1/ms)
	ah (1/ms) bh 	   (1/ms)
  hinf (1)
	df	(mV)
} 
STATE {
	m h
}
BREAKPOINT { 
	SOLVE states METHOD cnexp
	ina = gbar * m * m * m * h * ( v - ena ) 
	df = v - ena
} 
INITIAL { 
	settables( v )
	m = 0
	h  = hinf
} 
DERIVATIVE states { 
	settables( v ) 
	m' = am*( 1 - m ) - bm*m
	h' = (ah+bh)*(hinf - h)
}

UNITSOFF 

PROCEDURE settables(v(mV)) {

	TABLE am, bm, ah, bh, hinf  FROM -120 TO 40 WITH 641

	am = 0.182*(v+25)/(1-exp(-(v+25)/9))
  bm = -0.124*(v+25)/(1-exp((v+25)/9))
  ah = 0.024*(v+40)/(1-exp(-(v+40)/5))
  bh = -0.0091*(v+65)/(1-exp((v+65)/5))
  hinf = 1/(1+exp((v+55)/6.2))
}

UNITSON
