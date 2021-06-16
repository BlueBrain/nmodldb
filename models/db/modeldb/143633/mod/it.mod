TITLE McCormick and Huguenard low threshold calcium current

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX tcit
	NONSPECIFIC_CURRENT it
	RANGE pcatbar, m_inf, tau_m, h_inf, tau_h, qm, qh
	RANGE depth, cainf, taur
}

UNITS {
	(molar) = (1/liter)
	(mM)	= (millimolar)
	(mV) =	(millivolt)
	(mA) =	(milliamp)
	(um)	= (micron)

	FARADAY = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
}

PARAMETER {
	v		(mV)
	celsius		(degC)
	pcatbar	= 0.0000138	(cm/s)
	cao	= 2	(mM)
	qm	= 2.5
	qh	= 2.5
	depth	= .1	(um)
	taur	= 5	(ms)
	cainf	= 2.4e-4	(mM)
	vshift	= 3	(mV)
}

STATE {
	m h 
	cai (mM)
}

ASSIGNED {
	it	(mA/cm2)
	m_inf
	tau_m	(ms)
	h_inf
	tau_h	(ms)
	phi_m
	phi_h
	drive_channel	(mM/ms)
}

BREAKPOINT {
	SOLVE state METHOD euler
	it = pcatbar * m*m*h * ghk(v, cai, cao)
}

DERIVATIVE state {
	evaluate_fct(v)

	m' = (m_inf - m) / tau_m
	h' = (h_inf - h) / tau_h
	cai' = drive_channel + (cainf-cai)/taur
}


UNITSOFF
INITIAL {
	phi_m = qm ^ ((celsius-23.5)/10)
	phi_h = qh ^ ((celsius-23.5)/10)

	evaluate_fct(v)

	m = m_inf
	h = h_inf
	cai = cainf
}

PROCEDURE evaluate_fct(v(mV)) {
	
	drive_channel =  - (10000) * it / (2 * FARADAY * depth)
	if (drive_channel <= 0.) { drive_channel = 0. }	: cannot pump inward

:	m_inf = 1.0 / ( 1 + Exp(-(v+vshift+55)/6.2) ) : MB changed 57 ---> 52 to shift m_inf toward right
	m_inf = 1.0 / ( 1 + Exp(-(v+vshift+65)/6.2) )
: ORIG	m_inf = 1.0 / ( 1 + Exp(-(v+vshift+57)/6.2) )
:	h_inf = 1.0 / ( 1 + Exp((v+vshift+76)/4.0) ) : MB changed 81 --> 76 to shift h_inf toward right
	h_inf = 1.0 / ( 1 + Exp((v+vshift+81)/4.0) )

	tau_m = ( 0.612 + 1.0 / ( Exp(-(v+vshift+132)/16.7) + Exp((v+vshift+16.8)/18.2) ) ) / phi_m
	if( v < -80) {
		tau_h = Exp((v+vshift+467)/66.6) / phi_h
	} else {
		tau_h = 1.6147 * ( 28 + Exp(-(v+vshift+22)/10.5) ) / phi_h : 1.6147 entered by MB on 3/26/08 to make Tau_h curves match at -80.
	}
}

FUNCTION ghk(v(mV), ci(mM), co(mM)) (.001 coul/cm3) {
	LOCAL z, eci, eco
	z = (1e-3)*2*FARADAY*v/(R*(celsius+273.15))
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
		efun = z/(Exp(z) - 1)
	}
}

FUNCTION Exp(x) {
	if (x < -100) {
		:Exp = 0
	}else{
		if (x > 700) {
			Exp = exp(700)
		}else{
			Exp = exp(x)
		}
	}
}

UNITSON
