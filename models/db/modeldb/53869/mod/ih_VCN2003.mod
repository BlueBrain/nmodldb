TITLE Ih channels in VCN auditory neurons
: ih=ghmax*r*(v-eh) 
: based on Rothman and Manis (2003c)
: Modifications by Yi Zhou for an MSO model

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX ih_VCN2003
	USEION h READ eh WRITE ih VALENCE 1	
	RANGE ghbar 
	RANGE r_inf,tau_r,r_exp
	RANGE ih,gh
	
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	ghbar	= 0.002	(mho/cm2)  
	eh=-43		(mV) 
	celsius =22		(degC)
	dt              (ms)
	v               (mV)
	
}

STATE {
	r
}

ASSIGNED {
	gh	(mho/cm2)
	ih	(mA/cm2)
	r_inf
	tau_r
	r_exp
	tadj
}


BREAKPOINT {
	SOLVE states
	gh=ghbar *r
	ih  = gh*(v-eh)
}


PROCEDURE states() {	: this discretized form is more stable
	evaluate_fct(v)
	r = r + r_exp * (r_inf - r)
	VERBATIM
	return 0;
	ENDVERBATIM
}

UNITSOFF
INITIAL {
:
:  Q10 was assumed to be 3 for both currents
:
	tadj = 3.0 ^ ((celsius-22)/ 10 )
	evaluate_fct(v)
	r= r_inf
	}

PROCEDURE evaluate_fct(v(mV)) {
	
	tau_r = (100000/(237*exp((v+60)/12)+17*exp(-(v+60)/14))+25)/ tadj
	r_inf = 1/(1+exp((v+76)/7))
   
	r_exp = 1 - exp(-dt/tau_r)
	
}

UNITSON
