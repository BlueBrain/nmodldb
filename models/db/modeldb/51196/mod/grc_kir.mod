TITLE Cerebellum Granule Cell Model, Kir channel
COMMENT
Reference: E.D'Angelo, T.Nieus, A. Maffei, S. Armano, P. Rossi,
V. Taglietti, A. Fontana, G. Naldi "Theta-frequency bursting and 
resonance in cerebellar granule cells: experimental evidence and 
modeling of a slow K+-dependent mechanism", J. neurosci., 2001,
21,P. 759-770.
ENDCOMMENT
 
NEURON { 
	SUFFIX GrC_Kir 
	USEION k READ ek WRITE ik 
	RANGE gkbar, ik, g, alpha_d, beta_d 
	RANGE Aalpha_d, Kalpha_d, V0alpha_d
	RANGE Abeta_d, Kbeta_d, V0beta_d
	RANGE d_inf, tau_d 
} 
 
UNITS { 
	(mA) = (milliamp) 
	(mV) = (millivolt) 
} 
 
PARAMETER { 
	Aalpha_d = 0.13289 (/ms)
      Kalpha_d = -24.3902 (mV)
      V0alpha_d = -83.94 (mV)
	Abeta_d = 0.16994 (/ms)
	Kbeta_d = 35.714 (mV)
      V0beta_d = -83.94 (mV)
      gkbar = 0.0009 (mho/cm2) 
} 

STATE { 
	d 
} 

ASSIGNED { 
	ik (mA/cm2) 
	d_inf 
	tau_d (ms) 
	g (mho/cm2) 
	alpha_d (/ms) 
	beta_d (/ms) 
      ek (mV)
      celsius (degC) 
      v (mV) 
} 
 
INITIAL { 
	rate(v) 
	d = d_inf 
} 
 
BREAKPOINT { 
	SOLVE states METHOD derivimplicit
	g = gkbar*d   : primo ordine!!!
	ik = g*(v - ek) 
	alpha_d = alp_d(v) 
	beta_d = bet_d(v) 
} 
 
DERIVATIVE states { 
	rate(v) 
	d' =(d_inf - d)/tau_d 
} 
 
FUNCTION alp_d(v(mV))(/ms) { LOCAL Q10
	Q10 = 3^((celsius-20(degC))/10(degC))
if((v-V0alpha_d)/Kalpha_d >200){
alp_d = Q10*Aalpha_d*exp(200)
}else{
	alp_d = Q10*Aalpha_d*exp((v-V0alpha_d)/Kalpha_d) 
} 
} 
FUNCTION bet_d(v(mV))(/ms) { LOCAL Q10
	Q10 = 3^((celsius-20(degC))/10(degC))
if((v-V0beta_d)/Kbeta_d >200){
bet_d = Q10*Abeta_d*exp(200)
}else{
	bet_d = Q10*Abeta_d*exp((v-V0beta_d)/Kbeta_d) 
} 
 }
PROCEDURE rate(v (mV)) {LOCAL a_d, b_d 
	TABLE d_inf, tau_d  
	DEPEND Aalpha_d, Kalpha_d, V0alpha_d, 
	       Abeta_d, Kbeta_d, V0beta_d, celsius FROM -100 TO 100 WITH 200 
	a_d = alp_d(v)  
	b_d = bet_d(v) 
	tau_d = 1/(a_d + b_d) 
	d_inf = a_d/(a_d + b_d) 
} 

