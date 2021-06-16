TITLE Motor Axon Node channels

: 2/02  AXNODE.mod file
: McIntyre CC, Grill WM, Sherman DL, and Thakor NV. 2004. Cellular effects of deep brain : stimulation: model-based analysis of activation and inhibition. J Neurophysiol 91: 1457-1469.


NEURON {
	SUFFIX axnodenew2
	NONSPECIFIC_CURRENT ina
	NONSPECIFIC_CURRENT inap
	USEION k READ ek WRITE ik
	NONSPECIFIC_CURRENT il
	RANGE gnapbar, gnabar, gkbar, gl, ena, el
	RANGE mp_inf, m_inf, h_inf, s_inf
	RANGE tau_mp, tau_m, tau_h, tau_s
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gnapbar = 0.05	(mho/cm2)
	gnabar	= 3.0	(mho/cm2)
	gkbar   = 0.07 	(mho/cm2)
	gl	= 0.005 (mho/cm2)
	ena     = 45.0  (mV)
	el	= -70.0 (mV)
	hpChange  =  1
	mpChange  =  1
	mppChange  =  1
	spChange  =  1
	
}

STATE {
	mp m h s
}

ASSIGNED {
	ek (mV)
	v (mV)
	inap    (mA/cm2)
	ina	(mA/cm2)
	ik	(mA/cm2)
	il      (mA/cm2)
	mp_inf
	m_inf
	h_inf
	s_inf
	tau_mp (ms)
	tau_m (ms)
	tau_h (ms)
	tau_s (ms)
	q10_1
	q10_2
	q10_3
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	inap = gnapbar * mp*mp*mp * (v - ena)
	ina = gnabar * m*m*m*h * (v - ena)
	ik   = gkbar * s * (v - ek)
	il   = gl * (v - el)
}

DERIVATIVE states {   : exact Hodgkin-Huxley equations
       evaluate_fct(v)
	mp'= (mp_inf - mp) / tau_mp*mppChange
	m' = (m_inf - m) / tau_m*mpChange
	h' = (h_inf - h) / tau_h*hpChange
	s' = (s_inf - s) / tau_s*spChange
}

UNITSOFF

INITIAL {
	evaluate_fct(v)
	mp = mp_inf
	m = m_inf
	h = h_inf
	s = s_inf
}

PROCEDURE evaluate_fct(v(mV)) { LOCAL a,b

	a = vtrap1(v)
	b = vtrap2(v)
	tau_mp = 1 / (a + b)
	mp_inf = a / (a + b)

	a = vtrap6(v)
	b = vtrap7(v)
	tau_m = 1 / (a + b)
	m_inf = a / (a + b)

	a = vtrap8(v)
	b = 12.6 / (1 + Exp(-(v+21.8)/13.4))
	tau_h = 1 / (a + b)
	h_inf = a / (a + b)

	a = 0.3 / (Exp((v+43)/-5) + 1) 
	b = 0.03 / (Exp((v+80)/-1) + 1)
	tau_s = 1 / (a + b)
	s_inf = a / (a + b)
}

FUNCTION vtrap1(x) {
	if (fabs((x+17)/10.2) < 1e-6) {
		vtrap1 = 0.0353*10.2
	}else{
		vtrap1 = (0.0353*(x+17)) / (1 - Exp(-(x+17)/10.2))
	}
}

FUNCTION vtrap2(x) {
	if (fabs((x+24)/10) < 1e-6) {
		vtrap2 = 0.000883*10
	}else{
		vtrap2 = (0.000883*(-(x+24))) / (1 - Exp((x+24)/10))
	}
}

FUNCTION vtrap6(x) {
	if (fabs((x+11.4)/10.3) < 1e-6) {
		vtrap6 = 6.57*10.3
	}else{
		vtrap6 = (6.57*(x+11.4)) / (1 - Exp(-(x+11.4)/10.3))
	}
}

FUNCTION vtrap7(x) {
	if (fabs((x+15.7)/9.16) < 1e-6) {
		vtrap7 = 0.304*9.16
	}else{
		vtrap7 = (0.304*(-(x+15.7))) / (1 - Exp((x+15.7)/9.16))
	}
}

FUNCTION vtrap8(x) {
	if (fabs((x+104)/11) < 1e-6) {
		vtrap8 = 0.34*11
	}else{
		vtrap8 = (0.34*(-(x+104))) / (1 - Exp((x+104)/11)) 
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
