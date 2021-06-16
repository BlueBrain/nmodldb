TITLE Cerebellum Granule Cell Model

COMMENT
        Leakage

	Author: A. Fontana
	Last revised: 18.12.98
---
Adapted by Sungho Hong and Claus Lang
Computational Neuroscience Unit, Okinawa Institute of Science and Technology, Japan
Supervisor: Erik De Schutter

Correspondence: Sungho Hong (shhong@oist.jp)

September 16, 2017
ENDCOMMENT

NEURON {
	SUFFIX GRANULE_LKG1
	NONSPECIFIC_CURRENT il
	RANGE Q10_diff,g, fix_celsius
	RANGE el,ic, gbar_Q10, gbar
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	v (mV)
	gbar = 5.68e-5 (mho/cm2)
	Q10_diff	= 1.5
	el =  0 (mV)
    fix_celsius = 37 (degC)
}

ASSIGNED {
	il (mA/cm2)
	ic (mA/cm2)
	g (mho/cm2)
	gbar_Q10 (mho/cm2)
}

INITIAL {
  gbar_Q10 = gbar*(Q10_diff^((fix_celsius-30)/10))
  g = gbar_Q10
}
BREAKPOINT {
    il = g*(v - el)
    ic = il
}
