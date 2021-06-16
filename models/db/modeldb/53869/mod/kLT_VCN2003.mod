TITLE low threshold potassium channels in VCN auditory neurons 
: k_LT=glt*w^4*z*(v-Ek)
: based on Rothman and Manis (2003c)
: Modifications by Yi Zhou for an MSO model

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX kLT_VCN2003
	USEION k READ ek WRITE ik
	RANGE gkbar 
	RANGE w_inf,z_inf
	RANGE tau_w,tau_z
	RANGE w_exp,z_exp
	RANGE ik,gk
	
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gkbar	= 0.04	(mho/cm2)  
	ek=-70		(mV) 
	celsius =22		(degC)
	dt              (ms)
	v               (mV)
	
}

STATE {
	w z
}

ASSIGNED {
	gk(mho/cm2)
	ik(mA/cm2)
	w_inf
	z_inf
	tau_w
	tau_z
	w_exp
	z_exp
	tadj
}


BREAKPOINT {
	SOLVE states
	gk=gkbar *w^4*z
	ik  = gk*(v-ek)
}



PROCEDURE states() {	: this discretized form is more stable
	evaluate_fct(v)
	w = w + w_exp * (w_inf - w)
	z = z + z_exp * (z_inf - z)
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
	w= w_inf
	z= z_inf
	}

PROCEDURE evaluate_fct(v(mV)) {LOCAL h
	
	tau_w = (100/(6*exp((v+60)/6)+16*exp(-(v+60)/45))+1.5)/ tadj
	w_inf = 1/(1+exp(-(v+48)/6))^0.25
        
	tau_z = (1000/(exp((v+60)/20)+exp(-(v+60)/8))+50)/ tadj
	h=0.5
	z_inf = (1-h)/(1+exp((v+71)/10))+h
	
	w_exp = 1 - exp(-dt/tau_w)
	z_exp = 1 - exp(-dt/tau_z)
	
}

UNITSON
