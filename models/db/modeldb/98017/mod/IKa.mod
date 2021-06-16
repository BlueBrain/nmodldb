COMMENT
An adapting potassium current

Author: Fredrik Edin, 2003
Address: freedin@nada.kth.se

ENDCOMMENT

NEURON {
	SUFFIX IKa
	USEION k WRITE ik
	RANGE gk, ik, gkbar
	GLOBAL ainf, taua, binf, taub, ek
}

UNITS {
	(mV) = (millivolt)
	(mA) = (milliampere)
}

PARAMETER {
	gkbar 	= 0.001	(mho/cm2)	<0,1e9>
	ek 	= -80	(mV)
}

ASSIGNED {
	v	(mV)
	gk	(mho/cm2)
	ik	(mA/cm2)
	ainf
	taua	(ms)
	binf
	taub	(ms)	
}

STATE {
	a
	b
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	gk = gkbar * a^4 * b 
	ik = gk * ( v - ek )
}

DERIVATIVE state {
	rates( v )
	a' = ( ainf - a ) / taua
	b' = ( binf - b ) / taub
}


PROCEDURE rates( v (mV) ) {

	TABLE ainf, binf, taua, taub FROM -100 TO 100 WITH 200
	UNITSOFF
	ainf = 1 / ( 1 + exp( -(60 + v)/8.5 ) )
	taua = 0.37 + 1 / ( exp(-(v+238)/37.5) + exp((v+46)/5) ) 

	binf = 1 / ( 1 + exp((78 + v)/6 ) )	
	taub = 19 + 1 / ( exp(-(v+238)/37.5) + exp((v+46)/5) ) 
	UNITSON
} 
