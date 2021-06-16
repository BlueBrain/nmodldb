TITLE Motor Axon Node channels

: fast k in juxtaparanodal region model, based on:
:
: McIntyre CC, Grill WM, Sherman DL, and Thakor NV. 2004. Cellular effects of deep brain : stimulation: model-based analysis of activation and inhibition. J Neurophysiol 91: 1457-1469.


NEURON {
	SUFFIX fastK	
	USEION k READ ek WRITE ik
	NONSPECIFIC_CURRENT iflut
	RANGE gkfbar, gflut, eflut
	RANGE n_inf
	RANGE tau_n
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gkfbar = 0.04	(mho/cm2)
	gflut (mho/cm2)
	eflut=-70	(mV)
}

STATE {	n }

ASSIGNED {
	v (mV)
	ek (mV)
	ik      (mA/cm2)
	iflut   (mA/cm2)
	n_inf
	tau_n (ms)
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = gkfbar*n*n*n*n*(v - ek)
	iflut = gflut*(v-eflut)
}

DERIVATIVE states { 
        evaluate_fct(v)
	n'= (n_inf - n) / tau_n
}

UNITSOFF

INITIAL {
	evaluate_fct(v)
	n = n_inf
}

PROCEDURE evaluate_fct(v(mV)) { LOCAL a,b

	a = vtrap(v)
	b = vtrap0(v)
	tau_n = 1 / (a + b)
	n_inf = a / (a + b)

}

FUNCTION vtrap(x) {
	if (fabs((x+83.2)/1.1) < 1e-6) {
		vtrap = 0.0462*1.1
	}else{
		vtrap = (0.0462*(x+83.2)) / (1 - exp(-(x+83.2)/1.1))
	}
}

FUNCTION vtrap0(x) {
	if (fabs((x+66)/10.5) < 1e-6) {
		vtrap0 = 0.0824*10.5
	}else{
		vtrap0 = (0.0824*(-(x+66))) / (1 - exp((x+66)/10.5))
	}
}

UNITSON
