COMMENT
Sinusoidal afferent input from muscle spindle is modeled as a membrane mechanism.
ENDCOMMENT
					       
NEURON {
	SUFFIX IaSyn_Sinewave
	RANGE gmax, e, i, del, bias, g
	NONSPECIFIC_CURRENT i
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(S) = (siemens)
}

PARAMETER {
	gmax=0 	(S/cm2)	<0,1e9>
	e=0	(mV)
	
	del=0     (ms)
    tp=2000   (ms)
	bias=0 
	pie=3.14159
}

ASSIGNED {
	v (mV)
	i (mA/cm2)
	g (S/cm2)
}

BREAKPOINT {
	g = gmax * m(t)
	i = g*(v - e)
}

FUNCTION m(x) {
	at_time(del)

    if (t < del) {
       m=0   
    } else { 
        m = (sin(2*pie/tp*(t-tp/4-del)) + 1)/2
	}
}