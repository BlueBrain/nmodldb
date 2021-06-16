TITLE Cerebellum Granule Cell Model

COMMENT
        Calcium first order kinetics

	Author: A. Fontana
	Last revised: 12.12.98
---
Adapted by Sungho Hong and Claus Lang
Computational Neuroscience Unit, Okinawa Institute of Science and Technology, Japan
Supervisor: Erik De Schutter

Correspondence: Sungho Hong (shhong@oist.jp)

September 16, 2017
ENDCOMMENT

NEURON {
        SUFFIX GRANULE_CALC
        USEION ca READ ica, cao WRITE cai
	RANGE Q10_diff,beta_Q10, fix_celsius
        RANGE d, beta, cai0, ic
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
        ic             (mA/cm2)
        d = .2          (um)
        cao = 2.        (mM)
        cai0 = 1e-4     (mM)
	Q10_diff = 3
        beta = 1.5        (/ms)
        fix_celsius = 37 (degC)
}

ASSIGNED {
	beta_Q10 (mho/cm2)
  tau (ms)
}

STATE {
	cai (mM)
}

INITIAL {
	beta_Q10 = beta*(Q10_diff^((fix_celsius-30)/10))
        cai = cai0
        tau = (2*F*d)/(1e4)
}

BREAKPOINT {
    SOLVE conc METHOD derivimplicit
    ic = beta*(cai-cai0)
}

DERIVATIVE conc {
	cai' = -ica/tau - beta_Q10*(cai-cai0)
}
