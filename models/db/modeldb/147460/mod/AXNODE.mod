TITLE Motor Axon Node channels

: 2/02
: Cameron C. McIntyre
:
: Fast Na+, Persistant Na+, Slow K+, and Leakage currents 
: responsible for nodal action potential
: Iterative equations H-H notation rest = -80 mV
:
: This model is described in detail in:
:
: McIntyre CC, Richardson AG, and Grill WM. Modeling the excitability of
: mammalian nerve fibers: influence of afterpotentials on the recovery
: cycle. Journal of Neurophysiology 87:995-1006, 2002.

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX axnode
	NONSPECIFIC_CURRENT ina
	NONSPECIFIC_CURRENT inap
	NONSPECIFIC_CURRENT iks
	NONSPECIFIC_CURRENT ikf
	NONSPECIFIC_CURRENT il
	RANGE gnafbar, gnapbar, gksbar, gkfbar, gl, ena, ek, el
	RANGE m_inf, h_inf, p_inf, s_inf, n_inf
	RANGE tau_m, tau_h, tau_p, tau_s, tau_n
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {

	gnabar = 3.0	(mho/cm2)   :3.0
	gnapbar	= 0.01	(mho/cm2)
	gksbar   = 0.08 	(mho/cm2) :0.08
	gkfbar = 0.02 :0.02-0.04
	gl	= 0.007 (mho/cm2)
	ena= 50.0  (mV)  :50
	ek = -90.0 (mV)
	el	= -90.0 (mV)
	celsius		(degC)
	dt              (ms)
	v               (mV)
	
}

STATE {
	m h p s n
}

ASSIGNED {
	ina    (mA/cm2)
	inap	(mA/cm2)
	iks      (mA/cm2)
	ikf (mA/cm2)
	il      (mA/cm2)
	m_inf
	h_inf
	p_inf
	s_inf
	n_inf
	tau_m
	tau_h
	tau_p
	tau_s
	tau_n
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ina = gnabar * m*m*m*h * (v - ena)
	inap = gnapbar * p*p*p * (v - ena)
	iks   = gksbar * s* (v - ek)
	ikf = gkfbar*n*n*n*n*(v-ek)
	il   = gl * (v - el)
}

DERIVATIVE states {   : exact Hodgkin-Huxley equations
       evaluate_fct(v)
	m'= (m_inf - m) / tau_m
	h' = (h_inf - h) / tau_h
	p' = (p_inf - p) / tau_p
	s' = (s_inf - s) / tau_s
	n' = (n_inf - n) / tau_n
}

UNITSOFF

INITIAL {
:
:	Q10 adjustment
:



	evaluate_fct(v)
	m = m_inf
	h = h_inf
	p = p_inf
	s = s_inf
	n = n_inf
}

PROCEDURE evaluate_fct(v(mV)) { LOCAL a,b,v2

	a = vtrap1(v)
	b = vtrap2(v)
	tau_m = 1 / (a + b)
	m_inf = a / (a + b)

	a = vtrap3(v)
	b = vtrap4(v)
	tau_h = 1 / (a + b)
	h_inf = a / (a + b)
	
	a = vtrap5(v)
	b = vtrap6(v)
	tau_p = 1 / (a + b)
	p_inf = a / (a + b)
	
	a = vtrap7(v)
	b = vtrap8(v)
	tau_s = 1 / (a + b)
	s_inf = a / (a + b)
	
	a = vtrap9(v)
	b = vtrap10(v)
	tau_n = 1 / (a + b)
	n_inf = a / (a + b)

}


FUNCTION vtrap1(x) {
	if (fabs(-(x+20.4)/10.3) < 1e-6) {
		vtrap1 = 6.57*10.3
	}else{
		vtrap1 = (6.57*(x+20.4))/(1-Exp(-(x+20.4)/10.3))
	}
}

FUNCTION vtrap2(x) {
	if (fabs((x+25.7)/9.16) < 1e-6) {
		vtrap2 = 0.304*9.16
	}else{
		vtrap2 = (0.304*(-(x+25.7)))/(1-Exp((x+25.7)/9.16))  
	}
}

FUNCTION vtrap3(x) {
	if (fabs((x+114)/11) < 1e-6) {
		vtrap3 = 0.34*11
	}else{
		vtrap3 = (0.34*(-(x+114)))/(1-Exp((x+114)/11))
	}
}

FUNCTION vtrap4(x) {
		vtrap4 = 12.6/(1+Exp(-(x+31.8)/13.4))
	
}

FUNCTION vtrap5(x) {
	if (fabs(-(x+27)/10.2) < 1e-6) {
		vtrap5 = 0.0353*10.2
	}else{
		vtrap5 = (0.0353*(x+27)) / (1 - Exp(-(x+27)/10.2))
	}
}

FUNCTION vtrap6(x) {
	if (fabs((x+34)/10) < 1e-6) {
		vtrap6 = 0.000883*10 : Ted Carnevale minus sign bug fix
	}else{
		vtrap6 = (0.000883*(-(x+34))) / (1 - Exp((x+34)/10))
	}
}

FUNCTION vtrap7(x) {

		vtrap7 = 0.3 / (1 + Exp((x+53)/-5))
	
}

FUNCTION vtrap8(x) {

		vtrap8 = 0.03 / (1 + Exp((x+90)/-1))
	
}

FUNCTION vtrap9(x) {
	if (fabs(-(x+83.2)/1.1) < 1e-6) {
		vtrap9 = 0.0462*1.1 : Ted Carnevale minus sign bug fix
	}else{
		vtrap9 = (0.0462*(x+83.2))/(1-Exp(-(x+83.2)/1.1)) 
	}
}

FUNCTION vtrap10(x) {
	if (fabs((x+66)/10.5) < 1e-6) {
		vtrap10 = 0.0824*10.5
	}else{
		vtrap10 = (0.0824*(-(x+66)))/(1-Exp((x+66)/10.5))
	}
}



FUNCTION Exp(x) {
	if (x < -100) {
		Exp = 0
	}else{
		Exp = exp(x)
	}
}




UNITSON