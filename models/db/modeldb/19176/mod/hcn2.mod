COMMENT
Kinetic model of HCN2 channel gating from Wang et al 2002.

In this model channel opening is coupled to a change in the affinity of the cyclic nucleotide binding domain for cAMP which is manifest as a shift in the activation curve toward more positive potentials.  This model explains the slow activation kinetics of Ih associated with low concentrations of cAMP.

For further details email Matt Nolan at mfnolan@fido.cpmc.columbia.edu.

Reference

Wang J., Chen S., Nolan M.F. and Siegelbaum S.A. (2002). Activity-dependent regulation of HCN pacemaker channels by cyclicAMP: signalling through dynamic allosteric coupling. Neuron 36, 1-20.
ENDCOMMENT

NEURON {
	SUFFIX hcn2
	NONSPECIFIC_CURRENT i
	RANGE i, ehcn, g, gbar
	USEION a READ ai  VALENCE 0
}

UNITS {
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
}

PARAMETER {
	gbar = 1	(millimho/cm2)
	ehcn = -40 	(mV)
	a0 = 0.0015	(/ms)		: parameters for alpha and beta
	b0 = 0.02	(/ms)
	ah = -135.7    (mV)
	bh = -99.7    (mV)		
	ac = -0.155     (/mV)
	bc = 0.144     (/mV)
	aa0 = 0.0067	(/ms)		: parameters for alphaa and betaa
	ba0 = 0.014	(/ms)
	aah = -142.28   (mV)
	bah = -83.5   (mV)
	aac = -0.075    (/mV)
	bac = 0.144 (/mV)
	kon = 3085.7    (/mM-ms)		: cyclic AMP binding parameters
	koff = 0.000044857  (/ms)
	b  = 80
	bf = 8.94
	ai	(millimolar)      : concentration cyclic AMP
	gca = 1			: relative conductance of the bound state
	shift = 0	(mV)		: shift in voltage dependence
	q10v = 4                        : q10 value from Magee 1998
	q10a = 1.5			: estimated q10 for the cAMP binding reaction
	celsius         (degC)
}

ASSIGNED {
	v	(mV)
	g	(millimho/cm2)
	i	(milliamp/cm2)
	alpha	(/ms)
	beta    (/ms)
	alphaa	(/ms)
	betaa	(/ms)
}

STATE { c cac o cao }

INITIAL { SOLVE kin STEADYSTATE sparse }

BREAKPOINT {
	SOLVE kin METHOD sparse
	g = gbar*(o + cao*gca)
	i = g*(v-ehcn)*(1e-3)
}

KINETIC kin {
	LOCAL qa
	qa = q10a^((celsius-22 (degC))/10 (degC))
	rates(v)
	~ c <-> o       (alpha, beta)
	~ c <-> cac  (kon*qa*ai/bf,koff*qa*b/bf)
	~ o <-> cao      (kon*qa*ai, koff*qa)
	~ cac <-> cao      (alphaa, betaa)
	CONSERVE c + cac + o + cao = 1
}

PROCEDURE rates(v(mV)) {
	LOCAL qv
	qv = q10v^((celsius-22 (degC))/10 (degC))
	alpha = a0*qv / (1 + exp(-(v-ah-shift)*ac))
	beta = b0*qv / (1 + exp(-(v-bh-shift)*bc))
	alphaa = aa0*qv / (1 + exp(-(v-aah-shift)*aac))
	betaa = ba0*qv / (1 + exp(-(v-bah-shift)*bac))
}
