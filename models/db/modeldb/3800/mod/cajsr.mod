TITLE Calcium release from JSR
:from Courtemanche et al Am J Physiol 1998 275:H301


NEURON {	
SUFFIX Ca_jsr
	USEION cd WRITE icd VALENCE 2 : dummy ion to remove influence of icr on total membrane currents
	USEION cr READ icr, cri WRITE icr VALENCE 2 : ca-ion in JSR
	USEION ca READ cai, ica 
RANGE  Csqn, Kmcsqn, Vrel, krel
GLOBAL m_inf, h_inf, n_inf	
GLOBAL tau_m, tau_h, tau_n	
GLOBAL m_exp, h_exp, n_exp, Tauinactca	
GLOBAL dummy : prevent vectorization for use with CVODE

 }

UNITS {	
	(mA) = (milliamp)	
	(mV) = (millivolt)
	(mM) = (milli/liter)
        (mse-1) = (1/millisec)
}

PARAMETER {	
	krel = 30 (mse-1) <0,1e9>      
        Vrel = 96.48e-12 (cm3)
	Csqn = 10 (mM)
	Kmcsqn = 0.8 (mM)
	dt              (ms)	
	Tauinactca=1.91 (ms)
}

STATE {	m n h}

ASSIGNED {	v (mV)
	celsius (degC) : 37
	ica (mA/cm2)
        icr (mA/cm2)
	cai  (mM)
	cri  (mM)	
	icd  (mA/cm2)
	dummy
m_inf	h_inf	n_inf	tau_m	tau_h	tau_n	m_exp	h_exp	n_exp	tadj
}

BREAKPOINT {	SOLVE states METHOD derivimplicit 
	
	icr = krel*m*m*n*h*(cri - cai)
	icd = - icr
}
INITIAL {
evaluate_fct(icr,ica,v)	
m= m_inf	
n= n_inf
h= h_inf	
}
DERIVATIVE states {	
	evaluate_fct(icr,ica,v)
	m' = (m_inf - m)/tau_m	
	n' = (n_inf - n)/tau_n
	h' = (h_inf - h)/tau_h		
	


}

UNITSOFF


PROCEDURE evaluate_fct(icr,ica,v) { LOCAL a,b,q10	

q10 = 3.0 ^ ((celsius-37)/ 10 )
a = 1/(1 + expMa((Vrel*icr - ica/0.96480*100 - 3.4175e-1),(13.67e-4)))
b = 8/(q10*1(/ms))
tau_m = b
m_inf = a 
a = (1 - 1/(1 + expMb((Vrel*icr - ica/0.96480*100 - 6.835e-2),(13.67e-4)))) 
b = q10*1(/ms)*(Tauinactca + 2.09/(1 + expMa((Vrel*icr - ica/0.96480*100 - 3.4175e-1),(13.67e-4))))
tau_n = b
n_inf = a 
a = (1 - 1/(1 + exp(-(v - 40)/17)))
b = 6*(1 - exp(-(v - 7.9)/5))/((1 + 0.3*exp(-(v - 7.9)/5))*(v - 7.9))/(q10*1(/ms))
tau_h = b
h_inf = a 

}

FUNCTION expMa(x,y) {
	if (fabs(x/y) < 1e-6) {
		expMa = (1-x/y/2)
	}else{
		expMa = exp(-x/y)
	}
}

FUNCTION expMb(x,y) {
	if (fabs(x/y) < 1e-6) {
		expMb = (1-x/y/2)
	}else{
		expMb = exp(-x/y)
	}
}

UNITSON
