TITLE Thalamocortical neuron slow potassium channel

COMMENT

	Slow K+ current
	Implementation of Meijer et al., 2011
	Written by Xu Zhang, UConn, 2018

ENDCOMMENT

NEURON {
	SUFFIX tcslowK
	USEION k READ ek WRITE ik
	RANGE gbar, d_k, e1_k, e2_k, ekinf, ik
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
    gbar = 7e-4 (siemens/cm2)
	ek = -95 (mV)
}

ASSIGNED {
	v (mV)
	ik (mA/cm2)
	dinf
	taud (ms)
	
	ekinf
	tau1 (ms)
	tau2 (ms)
}

STATE {
	d_k
	e1_k
	e2_k
}

INITIAL {
	rates(v)
	d_k = dinf
	e1_k = ekinf
	e2_k = ekinf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = gbar * d_k * (0.4*e1_k+0.6*e2_k) * (v - ek)
}

DERIVATIVE states {
	rates(v)
	d_k' = (dinf-d_k)/taud 
	e1_k' = (ekinf-e1_k)/tau1
	e2_k' = (ekinf-e2_k)/tau2
}

PROCEDURE rates(v (mV)) {
	dinf = (1/(1+exp(-(v+43)/17)))^4
	taud = 2.5+0.253/(exp((v-81)/25.6)+exp(-(v+132)/18))
	
	ekinf = 1/(1+exp((v+58)/10.6))
	tau1  = 30.4+0.253/(exp((v-13.29)/200)+exp(-(v+130)/7.1))
	if (v<=-70) {
	tau2  = tau1
	} else {
	tau2 = 2260
	}
}