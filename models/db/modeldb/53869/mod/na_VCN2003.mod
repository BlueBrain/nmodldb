TITLE Na channels in VCN auditory neurons of guinea pig
 
: na=gna*m^3*h   
: based on Rothman and Manis 2003
: Modifications by Yi Zhou for an MSO model


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX na_VCN2003
	USEION na READ ena WRITE ina
	NONSPECIFIC_CURRENT il
	RANGE gnabar 
	RANGE m_inf,h_inf
	RANGE tau_m,tau_h
	RANGE m_exp,h_exp
	RANGE ina,gna
	RANGE gl,el
	
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gnabar	= 0.2	(mho/cm2) 
	ena=55		(mV) 
	gl= 4.0e-4	(mho/cm2)   
        el=-57		(mV)
	celsius=22 		(degC)
	dt              (ms)
	v               (mV)
	
}

STATE {
	m h
}

ASSIGNED {
	gna (mho/cm2)
	ina	(mA/cm2)
	il	(mA/cm2)
	m_inf
	h_inf
	tau_m
	tau_h
	m_exp
	h_exp
	tadj3
	
}


BREAKPOINT {
	SOLVE states
	gna=gnabar * m*m*h
	ina  = gna * (v - ena)
	il=gl*(v-el)
}



PROCEDURE states() {	: this discretized form is more stable
	evaluate_fct(v)
	m = m + m_exp * (m_inf - m)
	h = h + h_exp * (h_inf - h)
	VERBATIM
	return 0;
	ENDVERBATIM
}

UNITSOFF
INITIAL {
:
:  Q10 was assumed to be 3 for both currents
:
	tadj3 = 3.0 ^ ((celsius-22)/ 10 )
	
	evaluate_fct(v)
	m= m_inf
	h= h_inf
	}

PROCEDURE evaluate_fct(v(mV)) { 
	
	m_inf = 1 / (1+exp(-(v + 38) / 7))
    	h_inf = 1 / (1+exp((v + 65) / 6))

    	tau_m =  (10 / (5*exp((v+60) / 18) + 36*exp(-(v+60) / 25))) + 0.04
    	tau_h =  (100 / (7*exp((v+60) / 11) + 10*exp(-(v+60) / 25))) + 0.6
	
	tau_m=tau_m/tadj3	
	tau_h=tau_h/tadj3

	m_exp = 1 - exp(-dt/tau_m)
	h_exp = 1 - exp(-dt/tau_h)
	
}

UNITSON


