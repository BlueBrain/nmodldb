
COMMENT
a synaptic current with NMDA function conductance defined by
        i = g * (v - e)      i(nanoamps), g(micromhos);
        where
         g = 0 for t < onset and
         g = gmax * (exp(-(t - onset)/tau1)- exp(-(t - onset )/tau2))
                  / (1+etha*[Mg2+]*exp(-gamma*v))
          for t > onset
this has the property that the maximum value is gmax and occurs at
 t = delay + tau.
written by D.Rusakov, 1995
ENDCOMMENT
					       
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS NMDAsyn
	RANGE onset, tau1, tau2, Mgetha, gamma, gmax, e, i
	NONSPECIFIC_CURRENT i
}
UNITS {
	(molar) = (1/liter)
	(mM)	= (millimolar)
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
}

PARAMETER {
	onset=0 (ms)
	tau1=80 (ms)
	tau2=0.67 (ms)
	Mgetha=0.33
	gamma=0.06 (/mV)
	gmax=0 	(umho)
	e=0	(mV)
	v	(mV)
}

ASSIGNED { i (nA)  g (umho)}

BREAKPOINT {
	g = gmax*alpha( (t - onset)/tau1, (t - onset)/tau2)
	i = g*(v - e) / (1.+Mgetha*exp(-1.*gamma*v)) 
}

FUNCTION alpha(x1,x2) {
	if (x1 < 0 || x1> 10) { if (x2 < 0 || x2> 200) {
		alpha = 0
	}}else{
		alpha = exp(-x1)-exp(-x2)
	}
}

