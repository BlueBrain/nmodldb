NEURON {
	SUFFIX kv2_gp
	USEION k READ ek WRITE ik
	RANGE g, ik, an, bn, gbar
	GLOBAL a0, a1, ah, ac, b0, bc, Cq10
}

UNITS {
	(mV)	= (millivolt)
	(mA)	= (milliamp)
	(S)	= (siemens)
}

PARAMETER {
	gbar = 1	(S/cm2)
	celsius		(degC)
	Cq10 = 3


	a0 = 100.74	(1/ms)
	a1 = 0.612	(1/ms-mV)
	ah = -84.51	(mV)
	ac = -11.84	(mV)

	b0 = 0.0051	(1/ms)
	bc = 22.02	(mV)
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
	kf3	(1/ms)
	kb3	(1/ms)
	kf4	(1/ms)
	kb4	(1/ms)
}

STATE {
	c1
	c2
	c3
	c4
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
	~ c4 <-> c3     (kf1,kb1)
	~ c3 <-> c2     (kf2,kb2)
	~ c2 <-> c1     (kf3,kb3)
	~ c1 <-> o      (kf4,kb4)
	CONSERVE c4+c3+c2+c1+o=1
}

PROCEDURE rates(v(mV)) {LOCAL q10
	q10 = Cq10^((celsius-22 (degC))/10 (degC))
	an = q10*(a0 - a1 *v)/(exp((ah+v)/ac) - 1)
	bn = q10*(b0 )/(exp(v/bc))

	kf1 = 4*an
	kb1 = bn
	kf2 = 3*an
	kb2 = 2*bn
	kf3 = 2*an
	kb3 = 3*bn
	kf4 = an
	kb4 = 4*bn
}
