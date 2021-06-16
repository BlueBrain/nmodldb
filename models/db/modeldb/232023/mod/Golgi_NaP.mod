TITLE Cerebellum Golgi Cell Model

COMMENT
        pNa channel

	Author: E.D Angelo, T.Nieus, A. Fontana
	Last revised: 8.5.2000
---
Adapted by Sungho Hong and Claus Lang
Computational Neuroscience Unit, Okinawa Institute of Science and Technology, Japan
Supervision: Erik De Schutter

Correspondence: Sungho Hong (shhong@oist.jp)

September 16, 2017
ENDCOMMENT

NEURON {
	SUFFIX Golgi_NaP
	USEION na READ ena WRITE ina
	RANGE Q10_diff,Q10_channel,gbar_Q10, fix_celsius
	RANGE gbar, ina, g
	:RANGE Aalpha_m, Kalpha_m, V0alpha_m, alpha_m, beta_m
	:RANGE Abeta_m, Kbeta_m, V0beta_m
	:RANGE V0_minf, B_minf
	RANGE m, m_inf, tau_m, tcorr, ic
	:GLOBAL i
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gbar		= 0.00019 (mho/cm2)
	Aalpha_m 	= -0.91 (/mV-ms)
	Kalpha_m 	= -5 (mV)
	V0alpha_m 	= -40 (mV)
	Abeta_m 	= 0.62 (/mV-ms)
	Kbeta_m 	= 5 (mV)
	V0beta_m 	= -40 (mV)
	V0_minf 	= -43 (mV)
	B_minf 		= 5 (mV)
	v (mV)
	ena 	 (mV)
    fix_celsius = 37 (degC)
	Q10_diff	= 1.5
	Q10_channel	= 3
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
	tcorr	(1)
	gbar_Q10 (mho/cm2)
	ic
}

INITIAL {
	gbar_Q10 = gbar*(Q10_diff^((fix_celsius-23)/10))
	rate(v)
	m = m_inf
}

BREAKPOINT {
	SOLVE states METHOD derivimplicit
	g = gbar_Q10 * m
	ina = g*(v - ena)
	ic = ina
:	alpha_m = alp_m(v)
:	beta_m = bet_m(v)
}

DERIVATIVE states {
	rate(v)
	m' =(m_inf - m)/tau_m
}
FUNCTION alp_m(v(mV))(/ms) {
	tcorr = Q10_channel^((fix_celsius-30(degC))/10(degC))
	alp_m = tcorr * Aalpha_m*linoid(v-V0alpha_m, Kalpha_m)
}

FUNCTION bet_m(v(mV))(/ms) {
	tcorr = Q10_channel^((fix_celsius-30(degC))/10(degC))
	bet_m = tcorr * Abeta_m*linoid(v-V0beta_m, Kbeta_m)
}

PROCEDURE rate(v (mV)) {LOCAL a_m, b_m
	TABLE m_inf, tau_m
	DEPEND Aalpha_m, Kalpha_m, V0alpha_m,
	       Abeta_m, Kbeta_m, V0beta_m, fix_celsius FROM -100 TO 30 WITH 13000
	a_m = alp_m(v)
	b_m = bet_m(v)
:	m_inf = a_m/(a_m + b_m)
	m_inf = 1/(1+exp(-(v-V0_minf)/B_minf))
	tau_m = 5/(a_m + b_m)
}
FUNCTION linoid(x (mV),y (mV)) (mV) {
        if (fabs(x/y) < 1e-6) {
                linoid = y*(1 - x/y/2)
        }else{
                linoid = x/(exp(x/y) - 1)
        }
}
