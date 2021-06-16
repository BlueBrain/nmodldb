COMMENT
A slow potassium current

Author: Fredrik Edin, 2003
Address: freedin@nada.kth.se

ENDCOMMENT

NEURON {
	SUFFIX IKs
	USEION k READ ek WRITE ik
	RANGE gk, ik, ek, gkbar
	GLOBAL qinf, tauq, rinf, taur
}

UNITS {
	(mV) = (millivolt)
	(mA) = (milliampere)
}

PARAMETER {
	gkbar 	= 0.002	(mho/cm2)	<0,1e9>
}

ASSIGNED {
	v	(mV)
	ek	(mV)
	gk	(mho/cm2)
	ik	(mA/cm2)
	qinf
	tauq	(ms)
	rinf
	taur	(ms)	
}

STATE {
	r
	q
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	gk = gkbar * q * r 
	ik = gk * ( v - ek )
}

INITIAL {
	rates( v )
	q = qinf
	r = rinf
}	

DERIVATIVE state {
	rates( v )
	q' = ( qinf - q ) / tauq
	r' = ( rinf - r ) / taur
}

PROCEDURE rates( v (mV) ) {

        TABLE qinf, rinf, tauq, taur FROM -100 TO 100 WITH 200

	UNITSOFF
	qinf = 1 / ( 1 + exp( -(34 + v)/6.5 ) )
	tauq = 8 / ( exp(-(v+55)/30) + exp((v+55)/30) ) 

	rinf = 1 / ( 1 + exp((65 + v)/6.6 ) )	
	taur = 100 / ( 1 + exp(-(v+65)/6.8) ) + 100
	UNITSON
} 
