TITLE Cerebellum Granule Cell Model

COMMENT
        Gaba A leakage

	Author: A. Fontana
	Last revised: 18.2.99
---
Adapted by Sungho Hong and Claus Lang
Computational Neuroscience Unit, Okinawa Institute of Science and Technology, Japan
Supervisor: Erik De Schutter

Correspondence: Sungho Hong (shhong@oist.jp)

September 16, 2017
ENDCOMMENT

NEURON {
	SUFFIX GRANULE_LKG2
	NONSPECIFIC_CURRENT il
	RANGE Q10_diff,Q10_channel,gbar_Q10, fix_celsius
	RANGE egaba, g , ic, gbar
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	v (mV)
	gbar = 6e-5 (mho/cm2) : Increased of 200% for Jorntell
	egaba = -65 (mV)
	Q10_diff	= 1.5
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
    il = g*(v - egaba)
    ic =il
}
