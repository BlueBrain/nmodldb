: Six state HMM kinetic scheme for the WT Voltage-Gated Sodium Channel Nav1.7
: From the paper Kinetic Modeling of Nav1.7 Provides Insight Into Erythromelalgia-associated F1449V Mutation
: Gurkiewicz et al., J.Neurophysiol. (2011).
NEURON {
	SUFFIX wtnag3e
	USEION na READ ena WRITE ina
	RANGE g, gbar
}
UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
}

PARAMETER {
	gbar = 0			(pS/um2)
	a12 = 18			(/ms)
	a21 = 1.1			(/ms)
	a23 = 27			(/ms)
	a32 = 14			(/ms)
	a34 = 65			(/ms)
	a43 = 1				(/ms)
	a45 = 1.6			(/ms)
	a54 = 3e-5			(/ms)
	a56 = 51			(/ms)
	a65 = 3.8e8			(/ms)
	a35 = 0.22			(/ms)
	a53 = 6.3e-8		(/ms)
	a36 = 0.01			(/ms)
	a63 = 0.02			(/ms)
	a26 = 0.26			(/ms)
	a62 = 0.29			(/ms)
	z12 = 0.023			(/mV)
	z21 = 0.034			(/mV)
	z23 = 0.041			(/mV)
	z32 = 0.054			(/mV)
	z34 = 0.064			(/mV)
	z43 = 0.02			(/mV)
	z45 = 0.012			(/mV)
	z54 = 2e-5			(/mV)
	z56 = 0.035			(/mV)
	z65 = 0.22			(/mV)
	z35 = 0.057			(/mV)
	z53 = 0.039			(/mV)
	z36 = 0.095			(/mV)
	z63 = 0.184			(/mV)
	z26 = 1.4e-3		(/mV)
	z62 = 4.4e-3		(/mV)
}
ASSIGNED {
	v       (mV)
	ena     (mV)
	g       (pS/um2)
	ina     (mA/cm2)
	k12     (/ms)
	k21     (/ms)
	k23     (/ms)
	k32     (/ms)
	k34     (/ms)
	k43     (/ms)
	k45     (/ms)
	k54     (/ms)
	k56     (/ms)
	k65     (/ms)
	k35     (/ms)
	k53     (/ms)
	k36     (/ms)
	k63     (/ms)
	k26     (/ms)
	k62     (/ms)
}

STATE { c1 c2 c3 o I1 I2 }

BREAKPOINT {
	SOLVE states METHOD sparse
		g = gbar*o
		ina = (1e-4) * g*(v - ena)
}

INITIAL { SOLVE states STEADYSTATE sparse}

KINETIC states {
	rates(v)
        ~c1 <-> c2      (k12,k21)
        ~c2 <-> c3      (k23,k32)
        ~c3 <-> o       (k34,k43)
        ~o <-> I1       (k45,k54)
        ~I1 <-> I2      (k56,k65)
        ~c3 <-> I1      (k35,k53)
        ~c3 <-> I2      (k36,k63)
        ~c2 <-> I2      (k26,k62)
	CONSERVE c1+c2+c3+o+I1+I2=1
}

PROCEDURE rates(v(millivolt)) {

	k12 = a12*exp(z12*v)
	k21 = a21*exp(-z21*v)
	k23 = a23*exp(z23*v)
	k32 = a32*exp(-z32*v)
	k34 = a34*exp(z34*v)
	k43 = a43*exp(-z43*v)
	k45 = a45*exp(z45*v)
	k54 = a54*exp(-z54*v)
	k56 = a56*exp(z56*v)
	k65 = a65*exp(z65*v)
	k35 = a35*exp(z35*v)
	k53 = a53*exp(-z53*v)
	k36 = a36*exp(z36*v)
	k63 = a63*exp(z63*v)
	k26 = a26*exp(z26*v)
	k62 = a62*exp(-z62*v)
}