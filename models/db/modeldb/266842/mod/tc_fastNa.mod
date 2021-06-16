TITLE Thalamocortical neuron fast sodium channel

COMMENT

	Fast Na+ current
	Implementation of Meijer et al., 2011
	Written by Xu Zhang, UConn, 2018

ENDCOMMENT

NEURON {
	SUFFIX tcfastNa
	USEION na READ ena WRITE ina
	RANGE gbar, m_na, h_na, ina
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
    gbar = 3e-2 (siemens/cm2)
	ena = 45 (mV)
	
}

ASSIGNED {
	v (mV)
	ina (mA/cm2)
	minf
	taum (ms)
	alpham (1/ms)
	betam (1/ms)
	
	hinf
	tauh (ms)
	alphah (1/ms)
	betah (1/ms)
}

STATE {
	m_na
	h_na
}

INITIAL {
	rates(v)
	m_na = minf
	h_na = hinf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ina = gbar * (m_na^3) * h_na * (v - ena)
}

DERIVATIVE states {
	rates(v)
	m_na' = (minf-m_na)/taum 
	h_na' = (hinf-h_na)/tauh 
}

PROCEDURE rates(v (mV)) {
	alpham = 0.32*(v+55)/(1-exp(-(v+55)/4))
	betam = 0.28*(v+28)/(exp((v+28)/5)-1)
	minf = alpham/(alpham+betam) 
	taum = 1/(alpham + betam)
	
	alphah = 0.128*exp(-(v+51)/18)
	betah = 4/(1+exp(-(v+28)/5))
	hinf = alphah/(alphah+betah) 
	tauh = 1/(alphah + betah)
}