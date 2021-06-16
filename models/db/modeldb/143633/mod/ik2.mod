TITLE McCormick and Huguenard slow K channel

: M. Birdno changed proportion of h1 and h2 contributing to ik2
: in order to make it consistent with McCormick & Huguenard paper.

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX tcik2
	USEION k READ ek WRITE ik
	RANGE gk2bar, m_inf, tau_m, h1_inf, tau_h1, h2_inf, tau_h2, ek, i
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	v		(mV)
	celsius		(degC)
	dt		(ms)
:	ek	= -95	(mV)
	gk2bar	= 0.002	(mho/cm2)
	vshift	= 0	(mV)
}

STATE {
	m h1 h2
}

ASSIGNED {
	ek (mV)
	ik		(mA/cm2)
	i		(mA/cm2)
	m_inf
	tau_m		(ms)
	h1_inf
	tau_h1		(ms)
	h2_inf
	tau_h2		(ms)
	tadj
}

BREAKPOINT { 
	SOLVE states METHOD euler
	i = gk2bar * m * ((0.6 * h1)+(0.4 * h2)) * (v - ek) : Changed from 0.4*h1 and 0.6*h2 to match Huguenard paper
 	ik = i
}

DERIVATIVE states { 
	evaluate_fct(v)

	m'= (m_inf-m) / tau_m
	h1'= (h1_inf-h1) / tau_h1
	h2'= (h2_inf-h2) / tau_h2
}

UNITSOFF

INITIAL {
	tadj = 3^((celsius-23.5)/10)
	evaluate_fct(v)
	m = m_inf
     	h1 = h1_inf
	h2 = h2_inf
}

PROCEDURE evaluate_fct(v(mV)) {  LOCAL a,b
	tau_m = (1.0/(Exp((v+vshift-81)/25.6)+Exp((v+vshift+132)/-18))+9.9) / tadj
	m_inf = (1.0 / (1+Exp(-(v+vshift+43)/17)))^4


	tau_h1 = (1.0/(Exp((v+vshift-1329)/200)+Exp(-(v+vshift+129.7)/7.143))+120) / tadj
	h1_inf = 1.0/(1+Exp((v+vshift+58)/10.6))


	if (v<-70) {
		tau_h2 = tau_h1
		}
	else {
		tau_h2 = 8930 / tadj : Based on 8.9 s IK2b in Huguenard/McCormick paper
		}
	h2_inf = 1.0/(1+Exp((v+vshift+58)/10.6))
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