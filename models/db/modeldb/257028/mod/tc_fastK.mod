TITLE Thalamocortical neuron fast potassium channel

COMMENT

	Fast K+ current
	Implementation of Meijer et al., 2011
	Written by Xu Zhang, UConn, 2018

ENDCOMMENT

NEURON {
	SUFFIX tcfastK
	USEION k READ ek WRITE ik
	RANGE gbar, n_k, ik
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
    gbar = 3e-3 (siemens/cm2)
	ek = -95 (mV)
	
}

ASSIGNED {
	v (mV)
	ik (mA/cm2)
	ninf
	taun (ms)
	alphan (1/ms)
	betan (1/ms)
}

STATE {
	n_k
}

INITIAL {
	rates(v)
	n_k = ninf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = gbar * n_k^4 * (v - ek)
}

DERIVATIVE states {
	rates(v)
	n_k' = (ninf-n_k)/taun 
}

PROCEDURE rates(v (mV)) {
	alphan = 0.032*(v+63.8)/(1-exp(-(v+63.8)/5))
	betan = 0.5*exp(-(v+68.8)/40)
	ninf = alphan/(alphan+betan) 
	taun = 1/(alphan + betan)
}