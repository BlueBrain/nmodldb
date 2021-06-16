COMMENT
Described by Gruber et al. 2003.

Gruber, A.J., Solla, S.A., Surmeier, D.J., and Houk, J.C.
Modulation of striatal single units by expected reward: 
a spiny neuron model displaying dopamine-induced bistability.
J. Neurophysiol. 90:1095-1114, 2003.

The paper doesn't state the reversal potential for this current.
However, the figures are reproduced if e_leak equals e_k.
ENDCOMMENT

NEURON {
	SUFFIX leak
	NONSPECIFIC_CURRENT i
	RANGE g, i, iL
	RANGE e
}

UNITS {
	(mA) = (milliamp)
	(uA) = (microamp)
	(mV) = (millivolt)
	(mS) = (millimho)
}

PARAMETER {
	g = 0.008 (mS/cm2)	<0,1e9>
	e = -90	(mV)		: same as e_k in their model
}

ASSIGNED {
	v	(mV)
	iL	(uA/cm2)	: for consistency with their usage of uA/cm2
	i	(mA/cm2)
}

BREAKPOINT {
	iL = g*(v - e)
	i = (0.001)*iL
}
