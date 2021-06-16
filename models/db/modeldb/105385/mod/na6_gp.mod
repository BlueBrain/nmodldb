: NA6_GP.MOD
:
:				      	   b  -  isb01-isb02-isb03-isb04-isb05...-isb10
:				  		   |
: c1 - c2 - c3 - c4 - c5 - o  - iso01-iso02-iso03-iso04-iso05...-iso10
: |    |    |    |    |    |
: i1 - i2 - i3 - i4 - i5 - i6 - isi01-isi02-isi03-isi04-isi05...-isi10
:
: 01.24.2005
:
: Globus pallidus Nav1.6 channel

NEURON {
	SUFFIX na6_gp
	USEION na READ ena WRITE ina
	RANGE g, ina, gbar, m, bo, ob,Cn,Cf
	GLOBAL Con, Coff, n
	GLOBAL Oon0, hOon, cOon
	GLOBAL Ooff0, hOoff, cOoff
	GLOBAL a0, vha, vca
	GLOBAL b0, vhb, vcb
	GLOBAL g0
	GLOBAL d0
	GLOBAL aSo, aSi, aSb, bS, bS1, k
	GLOBAL SoR, SoH, SoC
	GLOBAL SiR, SiH, SiC
	GLOBAL SbR, SbH, SbC
	GLOBAL ob0, bl0, bslope, bvh
}

UNITS {
	(mV)	= (millivolt)
	(S)	= (siemens)
	(mA)	= (milliamp)
}

PARAMETER {
	gbar = 1		(S/cm2)

	a0 = 9			(1/ms)	: alpha
	vha  = -18		(mV)
	vca = 15		(mV)

	b0 = 5			(1/ms)	: beta
	vhb = -58		(mV)
	vcb = -60		(mV)

	g0 = 40			(1/ms)	: gamma

	d0 = 25			(1/ms)	: delta

	SoR = 0.0018	(1/ms)
	SoH	 = -25
	SoC  = 100

	SiR = 0.00023	(1/ms)
	SiH	 = -25
	SiC  = 100

	SbR = 0.00023	(1/ms)
	SbH	 = -25
	SbC  = 100

	bS1  = 0.0005	(1/ms)

	aSo = 3.3e-05	(1/ms)
	aSi = 3.3e-05	(1/ms)
	aSb = 3.3e-05	(1/ms)
	bS  = 4.4e-05		(1/ms)
	k = 4.0366

	Con = 0.025		(1/ms)
	Coff = 0.007	(1/ms)
	n = 2

	Oon0 = 0.77		(1/ms)
	hOon = -20		(mV)
	cOon = 30		(mV)

	Ooff0 = 0.001	(1/ms)
	hOoff = -20		(mV)
	cOoff = -50		(mV)

	ob0 = 1			(1/ms)	: "forward block rate"
	bl0 = 0.15		(1/ms)	: "reverse block rate"
	bslope = 25		(mV)	: "slope factor for block"
	bvh = 0			(mV)	: "half block voltage"

	Cq10 = 4
	celsius			(degC)
}

ASSIGNED {
	v	(mV)
	ena 	(mV)
	g	(S/cm2)
	ina	(mA/cm2)
	alpha  	(1/ms)
	beta	(1/ms)
	gamma   (1/ms)
	delta	(1/ms)
	Oon	(1/ms)
	Ooff	(1/ms)
	m
	ob	(1/ms)
	bo	(1/ms)
	Cn	(1/ms)
	Cf	(1/ms)
	q10
	aSo1
	aSi1
	aSb1
}

STATE {
	c1		: closed
	c2
	c3
	c4
	c5
	ct		: total closed
	o		: open
	i1		: fast inactivated
	i2
	i3
	i4
	i5
	i6
	ift		: total fast inactivated
	iso01	: slow inact from o
	iso02
	iso03
	iso04
	iso05
	iso06
	iso07
	iso08
	iso09
	iso10
	isot	: total slow inact from o
	isb01	: slow inact from b
	isb02
	isb03
	isb04
	isb05
	isb06
	isb07
	isb08
	isb09
	isb10
	isbt 	: total slow inact from b
	isi01	: slow inact from i6
	isi02
	isi03
	isi04
	isi05
	isi06
	isi07
	isi08
	isi09
	isi10
	isit 	: total slow inact from i6
	ist		: total slow inactivated
	it		: total inactivated
	bl		: (open) block
	avail	: fraction available
}

BREAKPOINT {
	SOLVE kin METHOD sparse
	g = gbar*o
	ina = g*(v-ena)
	ct = c1 + c2 + c3 + c4 + c5
	ift = i1 + i2 + i3 + i4 + i5 + i6
	isot = iso01 + iso02 + iso03 + iso04 + iso05 + iso06 + iso07 + iso08 + iso09 + iso10
	isit = isi01 + isi02 + isi03 + isi04 + isi05 + isi06 + isi07 + isi08 + isi09 + isi10
	isbt = isb01 + isb02 + isb03 + isb04 + isb05 + isb06 + isb07 + isb08 + isb09 + isb10
	ist = isot + isit + isbt
	it = ift + ist
	avail = o + ct
}

INITIAL {
	SOLVE kin STEADYSTATE sparse
}

KINETIC kin{
	rates(v)

	~ c1 <-> c2	(4*alpha, beta)
	~ c2 <-> c3	(3*alpha, 2*beta)
	~ c3 <-> c4	(2*alpha, 3*beta)
	~ c4 <-> c5	(alpha, 4*beta)
	~ c5 <-> o	(gamma, delta)
	~ o <-> iso01	(q10*aSo1, q10*bS1)
	~ iso01 <-> iso02	(q10*aSo*k, q10*bS*k)
	~ iso02 <-> iso03	(q10*aSo*k^2, q10*bS*k^2)
	~ iso03 <-> iso04	(q10*aSo*k^3, q10*bS*k^3)
	~ iso04 <-> iso05	(q10*aSo*k^4, q10*bS*k^4)
	~ iso05 <-> iso06	(q10*aSo*k^5, q10*bS*k^5)
	~ iso06 <-> iso07	(q10*aSo*k^6, q10*bS*k^6)
	~ iso07 <-> iso08	(q10*aSo*k^7, q10*bS*k^7)
	~ iso08 <-> iso09	(q10*aSo*k^8, q10*bS*k^8)
	~ iso09 <-> iso10	(q10*aSo*k^9, q10*bS*k^9)

	~ i1 <-> i2	(4*alpha*n, beta/n)
	~ i2 <-> i3	(3*alpha*n, 2*beta/n)
	~ i3 <-> i4	(2*alpha*n, 3*beta/n)
	~ i4 <-> i5	(alpha*n, 4*beta/n)
	~ i5 <-> i6	(gamma*m, delta/m)
	~ i6 <-> isi01	(q10*aSi1, q10*bS1)
	~ isi01 <-> isi02	(q10*aSi*k, q10*bS*k)
	~ isi02 <-> isi03	(q10*aSi*k^2, q10*bS*k^2)
	~ isi03 <-> isi04	(q10*aSi*k^3, q10*bS*k^3)
	~ isi04 <-> isi05	(q10*aSi*k^4, q10*bS*k^4)
	~ isi05 <-> isi06	(q10*aSi*k^5, q10*bS*k^5)
	~ isi06 <-> isi07	(q10*aSi*k^6, q10*bS*k^6)
	~ isi07 <-> isi08	(q10*aSi*k^7, q10*bS*k^7)
	~ isi08 <-> isi09	(q10*aSi*k^8, q10*bS*k^8)
	~ isi09 <-> isi10	(q10*aSi*k^9, q10*bS*k^9)

	~ c1 <-> i1	(Cn/n^4, Cf*n^4)
	~ c2 <-> i2	(Cn/n^3, Cf*n^3)
	~ c3 <-> i3	(Cn/n^2, Cf*n^2)
	~ c4 <-> i4	(Cn/n, Cf*n)
	~ c5 <-> i5	(Cn, Cf)
	~ o <-> i6	(Oon, Ooff)

	~ o <-> bl	(ob,bo)
	~ bl <-> isb01	(q10*aSb1, q10*bS1)
	~ isb01 <-> isb02	(q10*aSb*k, q10*bS*k)
	~ isb02 <-> isb03	(q10*aSb*k^2, q10*bS*k^2)
	~ isb03 <-> isb04	(q10*aSb*k^3, q10*bS*k^3)
	~ isb04 <-> isb05	(q10*aSb*k^4, q10*bS*k^4)
	~ isb05 <-> isb06	(q10*aSb*k^5, q10*bS*k^5)
	~ isb06 <-> isb07	(q10*aSb*k^6, q10*bS*k^6)
	~ isb07 <-> isb08	(q10*aSb*k^7, q10*bS*k^7)
	~ isb08 <-> isb09	(q10*aSb*k^8, q10*bS*k^8)
	~ isb09 <-> isb10	(q10*aSb*k^9, q10*bS*k^9)

	CONSERVE c1+c2+c3+c4+c5+i1+i2+i3+i4+i5+i6+iso01+iso02+iso03+iso04+iso05+iso06+iso07+iso08+iso09+iso10+isi01+isi02+isi03+isi04+isi05+isi06+isi07+isi08+isi09+isi10+isb01+isb02+isb03+isb04+isb05+isb06+isb07+isb08+isb09+isb10+o+bl=1


}

PROCEDURE rates(v(mV)) {
	q10 = Cq10^((celsius-23 (degC))/10 (degC))

	alpha = q10*a0*exp((v-vha)/vca)
	beta = q10*b0*exp((v-vhb)/vcb)
	gamma = q10*g0
	delta = q10*d0

	ob = q10*ob0
	bo = q10*bl0*exp(-(v-bvh)/bslope)

	Cn=q10*Con
	Cf=q10*Coff

	Oon = q10*Oon0*exp((v-hOon)/cOon)
	Ooff = q10*Ooff0*exp((v-hOoff)/cOoff)

	m = ((Oon/Ooff)/(Con/Coff))^(1/2)

	aSo1 = q10*SoR*exp((v-SoH)/SoC)
	aSi1 = q10*SiR*exp((v-SiH)/SiC)
	aSi1 = q10*SbR*exp((v-SbH)/SbC)
}
