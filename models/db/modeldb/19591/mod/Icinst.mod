: Calcium activated K channel.
: From Moczydlowski and Latorre (1983) J. Gen. Physiol. 82

UNITS {
	(molar) = (1/liter)
}

UNITS {
	(mV) =	(millivolt)
	(mA) =	(milliamp)
	(mM) =	(millimolar)
}


INDEPENDENT {t FROM 0 TO 1 WITH 100 (ms)}

NEURON {
	SUFFIX Ic
	USEION ca READ cai
	NONSPECIFIC_CURRENT i
	RANGE gbar
	GLOBAL oinf, tau
}

UNITS {
}

PARAMETER {
	v		(mV)
	gbar=.0873	(mho/cm2)	: Maximum Permeability
	cai       	(mM)
	e = -90		(mV)
	dt		(ms)
	

	f = 0.0851
	g = 0.077
	k1 = 1.5e-3	(mM)
	k2 = 1.5e-4	(mM)
	bbar = 1.5	(/ms)
	abar = 2.5	(/ms)
}
COMMENT
the preceding two numbers were switched on 8/19/92 in response to a bug
report by Bartlett Mel. In the paper the kinetic scheme is
C <-> CCa (K1)
CCa <-> OCa (beta2,alpha2)
OCa <-> OCa2 (K4)
In this model abar = beta2 and bbar = alpha2 and K4 comes from d2 and k2
I was forcing things into a nomenclature where alpha is the rate from
closed to open. Unfortunately I didn't switch the numbers.
ENDCOMMENT

ASSIGNED {
	i 		(mA/cm2)
	oinf
	tau		(ms)
}

STATE {	o }		: fraction of open channels

BREAKPOINT {
	SOLVE state
	i = gbar*o*(v - e)
}

LOCAL fac

:if state_cagk is called from hoc, garbage or segmentation violation will
:result because range variables won't have correct pointer.  This is because
: only BREAKPOINT sets up the correct pointers to range variables.
PROCEDURE state() {	: exact when v held constant; integrates over dt step
	rate(v, cai)
	o = o + fac*(oinf - o)
	VERBATIM
	return 0;
	ENDVERBATIM
}

INITIAL {
	rate(v, cai)
	o = oinf
}

FUNCTION alp(v (mV), ca (mM)) (1/ms) { :callable from hoc
	alp = abar/(1 + exp1(k1,f,v)/ca)
}

FUNCTION bet(v (mV), ca (mM)) (1/ms) { :callable from hoc
	bet = bbar/(1 + ca/exp1(k2,g,v))
}

FUNCTION exp1(k (mM), d, v (mV)) (mM) { :callable from hoc
	exp1 = k*exp(-d*v)
}

PROCEDURE rate(v (mV), ca (mM)) { :callable from hoc
	LOCAL a
	a = alp(v,ca)
	tau = 1/(a + bet(v, ca))
	oinf = a*tau
	fac = (1 - exp(-dt/tau))
}
