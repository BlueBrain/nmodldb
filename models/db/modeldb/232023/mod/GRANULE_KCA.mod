TITLE Cerebellum Granule Cell Model

COMMENT
        KCa channel

	Author: E.D'Angelo, T.Nieus, A. Fontana
	Last revised: 8.5.2000
---
Adapted by Sungho Hong and Claus Lang
Computational Neuroscience Unit, Okinawa Institute of Science and Technology, Japan
Supervisor: Erik De Schutter

Correspondence: Sungho Hong (shhong@oist.jp)

September 16, 2017
ENDCOMMENT

NEURON {
	SUFFIX GRANULE_KCA
	USEION k READ ek WRITE ik
	USEION ca READ cai
	RANGE Q10_diff,Q10_channel,gbar_Q10, fix_celsius
	RANGE gbar, ic, g, alpha_c, beta_c
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

	v (mV)
	cai (mM)
	Q10_diff	= 1.5
	Q10_channel	= 3
	gbar= 0.0045 (mho/cm2)
	ek = -84.69 (mV)
    fix_celsius = 37 (degC)
}

STATE {
	c
}

ASSIGNED {
	ic (mA/cm2)
	ik (mA/cm2)

	c_inf
	tau_c (ms)
	g (mho/cm2)
	alpha_c (/ms)
	beta_c (/ms)
	gbar_Q10 (mho/cm2)

  tcorr (1)

  bavc (mM)
  bbvc (mM)
}

INITIAL {
	gbar_Q10 = gbar*(Q10_diff^((fix_celsius-30)/10))
  tcorr = Q10_channel^((fix_celsius-30(degC))/10(degC))
	rate(v, cai)
	c = c_inf
}

BREAKPOINT {
	SOLVE states METHOD derivimplicit
	g = gbar_Q10*c
	ik = g*(v - ek)
	ic = ik
:  alp_c_bet_c(v, cai)
}

DERIVATIVE states {
	rate(v, cai)
	c' =(c_inf - c)/tau_c
}

PROCEDURE exprate(v(mv))(mM) {
  TABLE bavc, bbvc DEPEND Balpha_c, Kalpha_c, Bbeta_c, Kbeta_c FROM -100 TO 30 WITH 13000
  bavc = Balpha_c*exp(v/Kalpha_c)
  bbvc = Bbeta_c*exp(v/Kbeta_c)
}

PROCEDURE alp_c_bet_c(v(mV), cai(mM))(/ms) {
  exprate(v)
	alpha_c = tcorr*Aalpha_c/(1+(bavc/cai))
	beta_c = tcorr*Abeta_c/(1+cai/bbvc)
}

PROCEDURE rate(v (mV), cai(mM)) {
	alp_c_bet_c(v, cai)
	tau_c = 1/(alpha_c + beta_c)
	c_inf = alpha_c/(alpha_c + beta_c)
}
