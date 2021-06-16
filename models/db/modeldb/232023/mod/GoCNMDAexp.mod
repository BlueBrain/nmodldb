

NEURON {
	POINT_PROCESS GoCNMDAexp
	RANGE tau1, tau2, e, i,tau3,mg,alpha,beta,block
	NONSPECIFIC_CURRENT i

	RANGE g
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	tau1=  4  (ms) <1e-9,1e9> :50
	tau2 = 100 (ms) <1e-9,1e9> :165
	tau3 = 255 (ms):520
	e=0	(mV)
	mg		    = 1.2		(mM)
	alpha		= 0.062		(/mV)
	beta		= 3.57		(mM)
	block=1
}

ASSIGNED {
	v (mV)
	i (nA)
	g (uS)
	factor
}

STATE {
	A (uS)
	B (uS)
	C (uS)
}

INITIAL {
	
	
	LOCAL tp
     	if (tau1/tau2 > .9999) {
     		tau1 = .9999*tau2
     	}
     	
     	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
     	factor = -exp(-tp/tau1) + exp(-tp/tau2)
     	factor = 1/factor

	
	A = 0
	B = 0
	C = 0
	
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	block = 1/(1+(mg/beta)*exp(-alpha*v))

	g = block*((B+C) - A)
        i = g*(v - e)

}

DERIVATIVE state {
	A' = -A/tau1
	B' = -B/tau2
	C' = -C/tau3

}

NET_RECEIVE(weight (uS)) {
	A = A + weight*factor
	B = B + 0.75*weight*factor
	C = C + 0.25*weight*factor

}
