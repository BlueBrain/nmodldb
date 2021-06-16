TITLE Ca-dependent potassium current (C-current)

COMMENT
        *********************************************
        reference:      Yamada, Koch & Adams (1989) 
			Meth. in Neuronal Modeling, MIT press
        found in:       bullfrog sympathetic ganglion cells
        *********************************************
	Assembled for MyFirstNEURON by Arthur Houweling
Updated to use CVode - N-type VGCC dependent
ENDCOMMENT


NEURON {
	SUFFIX mykca
	USEION k READ ek WRITE ik
:	USEION ca READ cai
	USEION can READ cani
        RANGE gkbar, ik
	GLOBAL m_inf, tau_m
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mM) = (milli/liter)
}

PARAMETER {
	gkbar= 0.00345	(mho/cm2) 
}

ASSIGNED {
	v		(mV)
	celsius		(degC)
:	cai		(mM)
	cani     	(mM)
	ek		(mV)
	ik		(mA/cm2)
	tau_m		(ms)
	m_inf
	tadj
}

STATE {
	m
}

BREAKPOINT { 
	SOLVE states METHOD cnexp
	ik = gkbar * m * (v - ek)
}

DERIVATIVE states { 
	rates(v,cani)

       m'= (m_inf-m) / tau_m
}

INITIAL {
	tadj = 3^((celsius-23.5)/(10(degC)))
	rates(v,cani)
	m = m_inf
}

PROCEDURE rates( v(mV), cai(mM)) {  LOCAL a,b
	a = 250(/mM) * cai * exp(v/(24(mV))) 
	b = 0.1 * exp(-v/(24(mV)))
	tau_m = 1(ms)/(a+b) / tadj
	m_inf = a/(a+b)
}
