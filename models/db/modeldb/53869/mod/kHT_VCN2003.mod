TITLE high threshold potassium channels in VCN auditory neurons
  
: k_HT=ght*(rr*n^2+(1-rr)*p)*(v-Ek)
: based on Rothman and Manis (2003c)
:
: Modifications by Yi Zhou for an MSO model

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX kHT_VCN2003
	USEION k READ ek WRITE ik
	RANGE gkbar 
	RANGE n_inf,p_inf
	RANGE tau_n,tau_p
	RANGE n_exp,p_exp
	RANGE ik,gk
	
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gkbar	= 0.03	(mho/cm2)  
	ek=-70		(mV) 
	celsius =22		(degC)
	dt              (ms)
	v               (mV)
	
}

STATE {
	n p rr
}

ASSIGNED {
	gk(mho/cm2)
	ik	(mA/cm2)
	n_inf
	p_inf
	tau_n
	tau_p
	n_exp
	p_exp
	tadj
}


BREAKPOINT {
	SOLVE states
	gk=gkbar *(rr*n^2+(1-rr)*p)
	ik  = gk*(v-ek)
}



PROCEDURE states() {	: this discretized form is more stable
	evaluate_fct(v)
	n = n + n_exp * (n_inf - n)
	p = p + p_exp * (p_inf - p)
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
	n= n_inf
	p= p_inf
        rr=0.85
	}

PROCEDURE evaluate_fct(v(mV)) {
	
	tau_n = (100/(11*exp((v+60)/24)+21*exp(-(v+60)/23))+0.7)/ tadj
	n_inf = 1/(1+exp(-(v+15)/5))^2
        
	tau_p = (100/(4*exp((v+60)/32)+5*exp(-(v+60)/22))+5)/ tadj
	p_inf = 1/(1+exp(-(v+23)/6))
	
	n_exp = 1 - exp(-dt/tau_n)
	p_exp = 1 - exp(-dt/tau_p)
	
}

UNITSON
