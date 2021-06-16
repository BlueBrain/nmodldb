COMMENT
synaptic current with exponential rise and decay conductance defined by
        i = g * (v - e)      i(nanoamps), g(micromhos);
        where
         g = 0 for t < onset and
         g=amp*((1-exp(-(t-onset)/tau0))-(1-exp(-(t-onset)/tau1)))
          for t > onset
ENDCOMMENT
					       
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS synnm
	RANGE g, onset, tau0, tau1, gmax, e, i, mg
	NONSPECIFIC_CURRENT i
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
}

PARAMETER {
	onset=0  (ms)
	:tau1=86 (ms)		:Parameters from Destexhe et al (1994), through Abbott and Dayan
	:tau2=1.585 (ms)		:Parameters from Destexhe et al (1994), through Abbott and Dayan
	tau1=78 (ms)
	tau2=2 (ms)
	gmax=0 	 (umho)
	e=0	 (mV)
	v	 (mV)
	mg	= 1    (mM)		: external magnesium concentration
}

ASSIGNED { i (nA)  g (umho) }

LOCAL   a[2]
LOCAL   adjust

BREAKPOINT {
        g = cond(t)*mgblock(v)
	i = g*(v - e)
}

FUNCTION cond(x) {
	adjust=(tau1/tau2)^(tau2/(tau2-tau1))-(tau1-tau2)^(tau1/(tau2-tau1))
	if (x < onset) {
		cond = 0
	}else{
		a[0]=exp(-(x-onset)/tau1)
		a[1]=exp(-(x-onset)/tau2)
		cond = gmax*(a[0]-a[1])/adjust
	}
}

FUNCTION mgblock(v) {
		TABLE 
	DEPEND mg
	FROM -140 TO 80 WITH 1000

	: from Jahr & Stevens

	mgblock = 1 / (1 + exp(0.062 (/mV) * -v) * (mg / 3.57 (mM)))
}
