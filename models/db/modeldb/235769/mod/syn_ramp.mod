COMMENT
    Point process for generation of ascending and descending variation in synaptic conductance over time.
ENDCOMMENT 

NEURON {
	POINT_PROCESS RampSyn
	RANGE gmax, e, i, del, dur, pkamp, bias, g
	NONSPECIFIC_CURRENT i
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	gmax=0 	(uS)	<0,1e9>
	e=0	(mV)
	del=0   (ms)
    dur=20000   (ms)
    pkamp=1       
	bias=0 
}

ASSIGNED { v (mV) i (nA)  g (uS)}

BREAKPOINT {
	g = gmax * m(t)
	i = g*(v - e)
}

FUNCTION m(x) {
	at_time(del)
    at_time(del + dur)

    if (t < del) {
       m=0   
    } else { 
        if (t < del+dur/2) {
            m = (2*pkamp/dur)*t + bias
        } else{ 
            m = -(2*pkamp/dur)*(t-dur) + bias
		}
	}
}