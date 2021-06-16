NEURON {
	SUFFIX kv3_gp
	USEION k READ ek WRITE ik
	RANGE g, ik, an, bn, gbar
	GLOBAL a1, a2, a3:, a4
	GLOBAL b1, b2, b3
	GLOBAL Cq10
}

UNITS {
	(mV)	= (millivolt)
	(S)	= (siemens)
	(mA)	= (milliamp)
}

PARAMETER {
	gbar	= 1	(S/cm2)
	celsius		(degC)

	a1 = 0.05	(1/ms)
	a2 = 8		(mV)
	a3 = 8		(mV)

	b1 = 0.07	(1/ms)
	b2 = -12		(mV)
	b3 = -15		(mV)

	Cq10 = 3
}

ASSIGNED {
	v	(mV)
	ek	(mV)
	g	(S/cm2)
	ik	(mA/cm2)

	an	(1/ms)
	bn	(1/ms)
	kf1	(1/ms)
	kb1	(1/ms)
	kf2	(1/ms)
	kb2	(1/ms)
	q10
}

STATE {
	c
	o
}

BREAKPOINT {
	SOLVE kin METHOD sparse
	g = gbar*o
	ik = g*(v-ek)
}

INITIAL {
	SOLVE kin STEADYSTATE sparse
}

KINETIC kin{
	rates(v)
	~ c <-> o		(an, bn)
	CONSERVE	c+o=1
}

PROCEDURE rates(v(mV)) {
	q10 = Cq10^((celsius-22 (degC))/10 (degC))

	an = q10 * a1*exp((v-a2)/a3)
	bn = q10 * b1*exp((v-b2)/b3)

	kf1 = 2*an
	kb1 = bn
	kf2 = an
	kb2 = 2*bn
}