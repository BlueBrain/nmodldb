TITLE Thalamocortical neuron h channel

COMMENT

	h current
	Implementation of Meijer et al., 2011
	Written by Xu Zhang, UConn, 2018

ENDCOMMENT

NEURON {
	SUFFIX tch
	NONSPECIFIC_CURRENT ih
	RANGE gbar, m_h, ih
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
    gbar = 5e-4 (siemens/cm2)
	eh = -43 (mV)
	
}

ASSIGNED {
	v (mV)
	ih (mA/cm2)
	minf
	taum (ms)
	alpham (1/ms)
	betam (1/ms)
}

STATE {
	m_h
}

INITIAL {
	rates(v)
	m_h = minf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ih = gbar * (m_h^4) * (v - eh)
}

DERIVATIVE states {
	rates(v)
	m_h' = (minf-m_h)/taum 
}

PROCEDURE rates(v (mV)) {
	minf = 1/(1+exp((v+85)/5.5))
	taum = 1/(exp(-15.45-0.086*v)+exp(-1.17+0.0701*v))
}