: Six state HMM kinetic scheme for the F1449V Voltage-Gated Sodium Channel Nav1.7
: From the paper Kinetic Modeling of Nav1.7 Provides Insight Into Erythromelalgia-associated F1449V Mutation
: Gurkiewicz et al., J.Neurophysiol. (2011).

NEURON {
	SUFFIX fvnag3e
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
	a12 = 170			(/ms)
	a21 = 1.6			(/ms)
	a23 = 51			(/ms)
	a32 = 31			(/ms)
	a34 = 82			(/ms)
	a43 = 0.9			(/ms)
	a45 = 1.3			(/ms)
	a54 = 0.02			(/ms)
	a56 = 0.02			(/ms)
	a65 = 0.278			(/ms)
	a35 = 0.6			(/ms)
	a53 = 1e-4			(/ms)
	a36 = 4e-3			(/ms)
	a63 = 9.4e-6		(/ms)
	a26 = 0.7			(/ms)
	a62 = 1e-3			(/ms)
	z12 = 0.061			(/mV)
	z21 = 0.069			(/mV)
	z23 = 0.035			(/mV)
	z32 = 0.087			(/mV)
	z34 = 0.00733		(/mV)
	z43 = 0.029			(/mV)
	z45 = 0.025			(/mV)
	z54 = 9e-5			(/mV)
	z56 = 0.079			(/mV)
	z65 = 0.237			(/mV)
	z35 = 0.012			(/mV)
	z53 = 0.049			(/mV)
	z36 = 0.065			(/mV)
	z63 = 0.161			(/mV)
	z26 = 0.02			(/mV)
	z62 = 5.6e-3		(/mV)
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