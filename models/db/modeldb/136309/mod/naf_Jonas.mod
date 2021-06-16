TITLE Sodium transient current for Schmidt-Hieber et al. 2008, kinetics from Engel & Jonas, 2005.

COMMENT

	Implemented by Erin Munro
ENDCOMMENT

INDEPENDENT { t FROM 0 TO 1 WITH 1 (ms) }

UNITS { 
	(mV) = (millivolt) 
	(mA) = (milliamp) 
} 
NEURON { 
	SUFFIX nafJonas
	USEION na READ ena WRITE ina
	RANGE gbar, ina, m, h, df, am, bm, ah, bh
}
PARAMETER { 
	gbar = 0.0 	   (mho/cm2)
	v (mV) ena 		   (mV)  
} 
ASSIGNED { 
	ina 		   (mA/cm2) 
	am (1/ms) bm 	   (1/ms)
	ah (1/ms) bh 	   (1/ms)
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
	h  = ah/(ah+bh)
} 
DERIVATIVE states { 
	settables( v ) 
	m' = am*( 1 - m ) - bm*m
	h' = ah*( 1 - h ) - bh*h
}

UNITSOFF 

PROCEDURE settables(v(mV)) {

	TABLE am, bm, ah, bh  FROM -120 TO 40 WITH 641

	am = -93.8*(v - 117) / (exp( -(v - 117 ) / 17.7 ) - 1 )
  bm = 0.168*exp(-(v-12)/23.3)
  ah = 0.000354*exp(-(v-12)/18.7)
  bh = 6.63/(exp(-(v+5.68)/13.3)+1)
}

UNITSON
