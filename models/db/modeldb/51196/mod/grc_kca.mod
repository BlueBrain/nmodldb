TITLE Cerebellum Granule Cell Model, KCa channel

COMMENT
Reference: E.D'Angelo, T.Nieus, A. Maffei, S. Armano, P. Rossi,
V. Taglietti, A. Fontana, G. Naldi "Theta-frequency bursting and 
resonance in cerebellar granule cells: experimental evidence and 
modeling of a slow K+-dependent mechanism", J. neurosci., 2001,
21,P. 759-770.
ENDCOMMENT
 
NEURON { 
	SUFFIX GrC_KCa 
	USEION k READ ek WRITE ik 
	USEION ca READ cai
	RANGE gkbar, ik, ica, g, alpha_c, beta_c
	RANGE Aalpha_c, Balpha_c, Kalpha_c
	RANGE Abeta_c, Bbeta_c, Kbeta_c 
	RANGE c_inf, tau_c 
} 
 
UNITS { 
	(mA) = (milliamp) 
	(mV) = (millivolt) 
	(molar) = (1/liter)
	(mM) = (millimolar)
} 
 
PARAMETER { 
	Aalpha_c = 2.5 (/ms)
	Balpha_c = 1.5e-3 (mM)
	Kalpha_c =  -11.765 (mV)

	Abeta_c = 1.5 (/ms)
	Bbeta_c = 0.15e-3 (mM)
	Kbeta_c = -11.765 (mV)
      gkbar= 0.004 (mho/cm2)  
} 

STATE { 
	c 
} 

ASSIGNED { 
	ik (mA/cm2) 
	ica (mA/cm2)
      c_inf 
	tau_c (ms) 
	g (mho/cm2) 
	alpha_c (/ms) 
	beta_c (/ms) 
      ek (mV)
      celsius (degC) 
      v (mV) 
      cai (mM)
} 
 
INITIAL { 
	rate(v) 
	c = c_inf 
} 
 
BREAKPOINT { 
	SOLVE states METHOD derivimplicit 
	g = gkbar*c 
	ik = g*(v - ek) 
	alpha_c = alp_c(v) 
	beta_c = bet_c(v) 
} 
 
DERIVATIVE states { 
	rate(v) 
	c' =(c_inf - c)/tau_c 
} 
 
FUNCTION alp_c(v(mV))(/ms) { LOCAL Q10
	Q10 = 3^((celsius-30(degC))/10(degC))
if(v/Kalpha_c>200){ 
alp_c = Q10*Aalpha_c/(1+(Balpha_c*exp(200)/cai)) 
}else{
	alp_c = Q10*Aalpha_c/(1+(Balpha_c*exp(v/Kalpha_c)/cai)) 
} 
 }
FUNCTION bet_c(v(mV))(/ms) { LOCAL Q10
	Q10 = 3^((celsius-30(degC))/10(degC))
if(v/Kbeta_c>200){
bet_c = Q10*Abeta_c/(1+cai/(Bbeta_c*exp(200)))
}else{ 
	bet_c = Q10*Abeta_c/(1+cai/(Bbeta_c*exp(v/Kbeta_c))) 
} 
 }
PROCEDURE rate(v (mV)) {LOCAL a_c, b_c 
	a_c = alp_c(v)  
	b_c = bet_c(v) 
	tau_c = 1/(a_c + b_c) 
	c_inf = a_c/(a_c + b_c) 
} 

