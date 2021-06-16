TITLE somatic leak

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
   The leak conductance accounts for the
   somatic shunt (sharp electrode recording).
   The soma specific resistance is Rm = 2.5 kOhm.
   Another value must be used for dendritic membrane.

ENDCOMMENT

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
}

NEURON {
	SUFFIX WMPas2
	NONSPECIFIC_CURRENT i
	RANGE g, eleak
}

PARAMETER {
	v 		(mV)
	g = .0004	(mho/cm2)	<0,1e9>
	eleak = -65	(mV)
}

ASSIGNED { 
	i	(mA/cm2)
}

BREAKPOINT {
	i = g*(v - eleak)
}
