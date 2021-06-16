TITLE I-h channel from Kamilla Angelo, Michael London,Soren R. Christensen, and Michael Hausser

NEURON {
	SUFFIX h
	NONSPECIFIC_CURRENT i
        RANGE ghbar
        GLOBAL qinf,tauh,ratetau
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {
	v 		(mV)
      eh  =-40	(mV)        
	ghbar=.0001 	(mho/cm2)
	ratetau=1
}

STATE {
        q
}

ASSIGNED {
	i (mA/cm2)
      qinf      
      tauh
}

INITIAL {
	rate(v)
	q=qinf
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	i = q*ghbar*(v-eh)

}

DERIVATIVE states {  
	rate(v)
      q' =  (qinf - q)/tauh
}

PROCEDURE rate(v (mV)) {
	qinf=1/( 1+exp((90+v)/9.67) )
	:tauh=ratetau*1/0.00062*( exp((v+68)/-22) + exp((v+68)/7.14) )
      tauh=ratetau*1/(0.02*( exp((v+90)/-22) + exp((v+90)/22) ))

}














