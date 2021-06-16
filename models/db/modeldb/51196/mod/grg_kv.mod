TITLE Cerebellum Granule Cell Model, KDr channel

COMMENT   
Reference: E.D'Angelo, T.Nieus, A. Maffei, S. Armano, P. Rossi,
V. Taglietti, A. Fontana, G. Naldi "Theta-frequency bursting and 
resonance in cerebellar granule cells: experimental evidence and 
modeling of a slow K+-dependent mechanism", J. neurosci., 2001,
21,P. 759-770.
ENDCOMMENT

NEURON { 
	SUFFIX GrG_KV 
	USEION k READ ek WRITE ik 
	RANGE gkbar, ik, g, alpha_n, beta_n 
	RANGE Aalpha_n, Kalpha_n, V0alpha_n
	RANGE Abeta_n, Kbeta_n, V0beta_n
	RANGE n_inf, tau_n 
} 
 
UNITS { 
	(mA) = (milliamp) 
	(mV) = (millivolt) 
} 
 
PARAMETER { 
	Aalpha_n = -0.01 (/ms-mV)
	Kalpha_n = -10 (mV)
	V0alpha_n = -25 (mV)
	Abeta_n = 0.125 (/ms)
	Kbeta_n = -80 (mV)
	V0beta_n = -35 (mV)
	gkbar= 0.003 (mho/cm2) 
} 

STATE { 
	n 
} 

ASSIGNED { 
	ik (mA/cm2) 
	n_inf 
	tau_n (ms) 
	g (mho/cm2) 
	alpha_n (/ms) 
	beta_n (/ms) 
      ek(mV)
      celsius (degC) 
      v (mV) 
} 
 
INITIAL { 
	rate(v) 
	n = n_inf 
} 
 
BREAKPOINT { 
	SOLVE states METHOD derivimplicit 
	g = gkbar*n*n*n*n 
	ik = g*(v - ek) 
	alpha_n = alp_n(v) 
	beta_n = bet_n(v) 
} 
 
DERIVATIVE states { 
	rate(v) 
	n' =(n_inf - n)/tau_n 
} 
 
FUNCTION alp_n(v(mV))(/ms) { LOCAL Q10
	Q10 = 3^((celsius-6.3(degC))/10(degC)) 
	alp_n = Q10*Aalpha_n*linoid(v-V0alpha_n, Kalpha_n)
} 
 
FUNCTION bet_n(v(mV))(/ms) { LOCAL Q10
	Q10 = 3^((celsius-6.3(degC))/10(degC)) 
if((v-V0beta_n)/Kbeta_n>200){
bet_n = Q10*Abeta_n*exp(200)
}else{
	bet_n = Q10*Abeta_n*exp((v-V0beta_n)/Kbeta_n) 
} 
} 
PROCEDURE rate(v (mV)) {LOCAL a_n, b_n 
	TABLE n_inf, tau_n 
	DEPEND Aalpha_n, Kalpha_n, V0alpha_n, 
               Abeta_n, Kbeta_n, V0beta_n, celsius FROM -100 TO 100 WITH 200 
	a_n = alp_n(v)  
	b_n = bet_n(v) 
	tau_n = 1/(a_n + b_n) 
	n_inf = a_n/(a_n + b_n) 
} 

FUNCTION linoid(x (mV),y (mV)) (mV) {
        if (fabs(x/y) < 1e-6) {
                linoid = y*(1 - x/y/2)
        }else{
if(x/y>200){
linoid = x/(exp(200) - 1)
}else{
                linoid = x/(exp(x/y) - 1)
        }
}
}

