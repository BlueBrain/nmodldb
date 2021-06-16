TITLE Thalamocortical neuron T-type calcium channel

COMMENT

	T-type calcium channel
	Implementation of Meijer et al., 2011
	Written by Xu Zhang, UConn, 2018

ENDCOMMENT

NEURON {
	SUFFIX tcCaT
	USEION ca READ cai, cao WRITE ica
	RANGE pca, m_t, minf, taum, h_t, hinf, tauh, cai, pump, ica
}

UNITS {
	(molar) = (1/liter)
	(mV) =	(millivolt)
	(mA) =	(milliamp)
	(mM) =	(millimolar)
	(um)	= (micron)
	FARADAY = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
}

PARAMETER {
	v		(mV)
	pca	= 1e-4	(cm/s)
	cao	= 2	(mM)
	dpt  = 1.0 (um)
	tau_ca = 5.0 (ms)
	kca = 5.1821e-5
}

STATE {
	m_t
	h_t
}

ASSIGNED {
	ica	(mA/cm2)
	cai (mM)
	minf
	taum	(ms)
	hinf
	tauh	(ms)
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	
	ica = pca*(m_t^2)*h_t*ghk(v,cai,cao)
}

DERIVATIVE states {
	rates(v)
	m_t' = (minf-m_t)/taum 
	h_t' = (hinf-h_t)/tauh
}

INITIAL {
	m_t = 0
	h_t = 0
	:cai = cabuf
}

PROCEDURE rates(v (mV)) {
	minf = 1/(1+exp(-(v+60)/6.2))
	taum = 0.204+0.333/(exp(-(v+136)/16.7)+exp((v+19.8)/18.2))
	
	hinf = 1/(1+exp((v+84)/4))
	if (v>=-81) {
	tauh  = 9.33+0.333*exp(-(v+25)/10.5)
	} else {
	tauh = 0.333*exp((v+470)/66.6)
	}
}

FUNCTION ghk(v(mV), ci(mM), co(mM)) (.001 coul/cm3) {
	LOCAL z, eci, eco
	z = (1e-3)*2*FARADAY*v/(R*(36+273.15))
	eco = co*efun(z)
	eci = ci*efun(-z)
	:high cao charge moves inward
	:negative potential charge moves inward
	ghk = (.001)*2*FARADAY*(eci - eco)
}

FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}