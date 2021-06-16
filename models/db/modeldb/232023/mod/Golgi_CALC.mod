TITLE Cerebellum Golgi Cell Model

COMMENT
        Calcium first order kinetics

	Author: A. Fontana
	Last revised: 12.12.98
---
Adapted by Sungho Hong and Claus Lang
Computational Neuroscience Unit, Okinawa Institute of Science and Technology, Japan
Supervision: Erik De Schutter

Correspondence: Sungho Hong (shhong@oist.jp)

September 16, 2017
ENDCOMMENT

NEURON {
        SUFFIX Golgi_CALC
        USEION ca READ ica, cao WRITE cai
	RANGE Q10_diff,beta_Q10, ic, fix_celsius
        RANGE d, beta, cai0, ca_pump_i
}

UNITS {
        (mV)    = (millivolt)
        (mA)    = (milliamp)
	(um)    = (micron)
	(molar) = (1/liter)
        (mM)    = (millimolar)
   	F      = (faraday) (coulomb)
}

PARAMETER {
        ica             (mA/cm2)
        fix_celsius = 37 (degC)
        d = .2          (um)
        cao = 2.        (mM)
        cai0 = 1e-4     (mM)
        beta = 1.3        (/ms)
	Q10_diff = 1.7
}

ASSIGNED {
	beta_Q10 (mho/cm2)
	ca_pump_i	(mA)
	ic
  tau (ms)
}
STATE {
	cai (mM)
}

INITIAL {
	beta_Q10 = beta*(Q10_diff^((fix_celsius-23)/10))
        cai = cai0
        tau = (2*F*d)/(1e4)
}

BREAKPOINT {
       SOLVE conc METHOD derivimplicit
}

DERIVATIVE conc {
	:  outward ionic current with valence 2+
	ca_pump_i = 2*beta_Q10*(cai-cai0)
	ic = ca_pump_i
	:  total outward Ca current
	cai' =  -ica/tau - beta_Q10*(cai-cai0)
}
