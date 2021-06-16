TITLE Cerebellum Granule Cell Model, pNa channel

COMMENT
Reference: E.D'Angelo, T.Nieus, A. Maffei, S. Armano, P. Rossi,
V. Taglietti, A. Fontana, G. Naldi "Theta-frequency bursting and 
resonance in cerebellar granule cells: experimental evidence and 
modeling of a slow K+-dependent mechanism", J. neurosci., 2001,
21,P. 759-770.
ENDCOMMENT
 
NEURON { 
	SUFFIX GrC_pNa 
	USEION na READ ena WRITE ina 
	RANGE gnabar, ina, g, alpha_m, beta_m
	RANGE Aalpha_m, Kalpha_m, V0alpha_m
	RANGE Abeta_m, Kbeta_m, V0beta_m
	RANGE V0_minf, B_minf
	RANGE m_inf, tau_m
} 
 
UNITS { 
	(mA) = (milliamp) 
	(mV) = (millivolt) 
} 
 
PARAMETER { 
	Aalpha_m = -0.091 (/mV-ms) 
	Kalpha_m = -5 (mV)
	V0alpha_m = -42 (mV)
	Abeta_m = 0.062 (/mV-ms)  
	Kbeta_m = 5 (mV)
	V0beta_m = -42 (mV)
	V0_minf = -42 (mV)
	B_minf = 5 (mV)
	gnabar= 2e-5 (mho/cm2)  
} 

STATE { 
	m 
} 

ASSIGNED { 
	ina (mA/cm2) 
	m_inf 
	tau_m (ms) 
	g (mho/cm2) 
	alpha_m (/ms)
	beta_m (/ms)
      ena (mV)
      celsius (degC) 
      v (mV) 
} 
 
INITIAL { 
	rate(v) 
	m = m_inf 
} 
 
BREAKPOINT { 
	SOLVE states METHOD derivimplicit 
	g = gnabar*m 
	ina = g*(v - ena) 
	alpha_m = alp_m(v)
	beta_m = bet_m(v)
} 
 
DERIVATIVE states { 
	rate(v) 
	m' =(m_inf - m)/tau_m 
} 

FUNCTION alp_m(v(mV))(/ms) { LOCAL Q10
	Q10 = 3^((celsius-30(degC))/10(degC))
	alp_m = Q10 * Aalpha_m*linoid(v-V0alpha_m, Kalpha_m) 
} 
 
FUNCTION bet_m(v(mV))(/ms) { LOCAL Q10
	Q10 = 3^((celsius-30(degC))/10(degC))
	bet_m = Q10 * Abeta_m*linoid(v-V0beta_m, Kbeta_m) 
} 
 
PROCEDURE rate(v (mV)) {LOCAL a_m, b_m 
	TABLE m_inf, tau_m 
	DEPEND Aalpha_m, Kalpha_m, V0alpha_m, 
	       Abeta_m, Kbeta_m, V0beta_m, celsius FROM -100 TO 100 WITH 200 
	a_m = alp_m(v)  
	b_m = bet_m(v)  
if((-(v-V0_minf)/B_minf)>200){
m_inf = 1/(1+exp(200))
}else{
	m_inf = 1/(1+exp(-(v-V0_minf)/B_minf))
}
	tau_m = 5/(a_m + b_m) 
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


