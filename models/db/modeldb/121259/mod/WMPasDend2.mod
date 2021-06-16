TITLE dendritic leak 

COMMENT

------------------------------------------
Ragnhild Halvorsrud June 1999
Dep. of Physiology, University of Oslo
------------------------------------------

Reference:
Lyle Borg-Graham:  'Interpretations of Data and Mechanisms for HPC-Models'
in "Cerebral Cortex" vol. 12: Cortical Models.
Eds: Jones, E. G. and Ulinski, P. S.

      Q10: Temperature sensitivity
    T_ref: Reference temperature
   T_corr: Correction factor when using other temperatures than T_ref

NB NB 
   The leak conductance represents the dendrite membrane property.
   The dendrite specific resistance is Rm = 40 kOhm.
   Another value must be used for somatic membrane.

ENDCOMMENT


UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
}

NEURON {
	SUFFIX WMPasDend2
	NONSPECIFIC_CURRENT i
	RANGE g, eleak
}

PARAMETER {
	v		(mV)
	g = 2.5e-5	(mho/cm2)	<0,1e9>
	eleak = -65	(mV)
}

ASSIGNED { 
	i	(mA/cm2)
}

BREAKPOINT {
	i = g*(v - eleak)
}
