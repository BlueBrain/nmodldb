
TITLE Glutmate Synaptic current

COMMENT

This file is based on exp2nmdar.mod of Baker et al. J Comp Nsc 2011.
It was modified to provide Glutamatergic synaptic current that is comprised of 
1) an alpha function model of AMPA-R
2) a bi-exponential model of an NMDA-R receptor synapse with a Jahr&Stevens Mg++ voltage dependency.
In addition, SHL newly implemented this file with the model of short-term facilitation and depression.

For Depression
Let D(n)_ = RRP just before nth AP
Let _D(n) =  RRP just after nth AP
_Dn = Dn_ (1- Pn) = Dn_ (1- Pb*F)
D(n+1)_ = _Dn + (1 - _Dn)(1-exp(-dt/tauD))

For Facilitation
Let Pb = basal Pr.
Let Af = inc in Pr just after an AP.
P(n+1) = Pb + (dPn + Af) exp(-dt/tauf)
Let Fn = Pn/Pb.
Since dPn/Pb = Fn - 1, F(n+1) = 1 + (Fn - 1 + Af/Pb) exp(-dt/tauf)
Let Af/Pb = f.
F(n+1) = 1 + (Fn + f - 1) exp(-dt/tauf)

******************************************************************************************
Here, is the comments of Baker et al. (2011).

This provides a simple dual-exponential model of an NMDA receptor synapse with a Jahr&Stevens Mg++ voltage dependency.

Changes were made by John Baker to the standard exp2syn.mod file so that the voltage dependency is addressed. 
The mgblock code was borrowed from a model by A. Destexhe.

The NMDA receptor is temperature sensitive. Any necesary adjustment to the time constants should be done by
setting tau1 and tau2 via hoc.

Default values are more or less taken from Dalby and Mody,J Neurophysiol 90: 786-797, 2003. No strong claims for physiological accuracy are made here.

--- (and now back to the original exp2syn comments) --------------------------------------

Two state kinetic scheme synapse described by rise time tau1,and decay time constant tau2. The normalized peak condunductance is 1.
Decay time MUST be greater than rise time.

The solution of A->G->bath with rate constants 1/tau1 and 1/tau2 is
 A = a*exp(-t/tau1) and
 G = a*tau2/(tau2-tau1)*(-exp(-t/tau1) + exp(-t/tau2)), where tau1 < tau2

If tau2-tau1 -> 0 then we have a alphasynapse.
and if tau1 -> 0 then we have just single exponential decay.

The factor is evaluated in the initial block such that an event of weight 1 generates a peak conductance of 1.

Because the solution is a sum of exponentials, the coupled equations can be solved as a pair of independent equations
by the more efficient cnexp method. 
******************************************************************************************

ENDCOMMENT

: Declare public variables
NEURON {
	POINT_PROCESS GluSyn
	RANGE ntar, e, i, mg, mgshift, tau1, tau2, tau3, tauD, tauF, f, Pb, Gnmda, Gampa
	NONSPECIFIC_CURRENT i, inmda, iampa
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
	(mM) = (milli/liter)
}

PARAMETER {
	tau1 = 0.5 (ms)			  : tau of alpha synapse for AMPA-R
	tau2 = 4 (ms) <1e-9,1e9>  : or 4 ms (Baker), rise time of Inmda
	tau3 = 42 (ms) <1e-9,1e9> : 42 (Baker) or 70 ms (Larkum), decay time of Inmda
	e = 0.0	(mV)
	mg = 1 (mM) 			  : external magnesium concentration
	sh = 0 (mV)
	ntar = .3	(1) < 0, 1 >  : NMDA to AMPA ratio
	
	f = 2 (1) < 0, 1e9 >    : Pr inc factor just after an AP
	tauF = 100 (ms) < 1e-9, 1e9 > : decay tau of f
	tauD = 500 (ms) < 1e-9, 1e9 > : RRP recovery tau
	Pb = 0.3 (1) < 0, 1 >	
}

ASSIGNED {
	v (mV)
	i (nA)
	inmda (nA)
	iampa (nA)
	Gnmda (uS)
	factor
}

STATE {
	Anmda (uS)
	Bnmda (uS)
	Aampa (uS)
	Gampa (uS)
}

INITIAL {
	LOCAL tp
	if (tau2/tau3 > .9999) {
		tau2 = .9999*tau3
	}
	Aampa = 0
	Gampa = 0
	Anmda = 0
	Bnmda = 0
	tp = (tau2*tau3)/(tau3 - tau2) * log(tau3/tau2)
	factor = -exp(-tp/tau2) + exp(-tp/tau3)
	factor = 1/factor
	mgblock(v)
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	Gnmda = Bnmda - Anmda
	inmda = Gnmda*(v - e)*mgblock(v)
	iampa = Gampa *(v - e)
	i = inmda + iampa
}

DERIVATIVE state {
	Aampa' = - Aampa/tau1
	Gampa' = Aampa/tau1 - Gampa/tau1	: Aampa -> Gampa -> disappear with rate const of 1/tau1.
	
	Anmda' = -Anmda/tau2
	Bnmda' = -Bnmda/tau3
}

NET_RECEIVE(weight (uS), F, D, tsyn (ms)) {
	INITIAL {
		: This header appears once per a NetStim event (stream)
		F = 1
		D = 1
		tsyn = t
    }
	F = 1 + (F-1)*exp(-(t - tsyn)/tauF)
	D = 1 - (1-D)*exp(-(t - tsyn)/tauD)
	tsyn = t 	: store previous stim time
	
	Aampa = Aampa + weight*exp(1)*F*D*Pb
	Anmda = Anmda + ntar*weight*factor*F*D*Pb
	Bnmda = Bnmda + ntar*weight*factor*F*D*Pb

	F = F + f  : Pr just after an AP
	D = D*(1 - Pb*F)	: RRP just after an AP
}

: The following is borrowed from Destexhe NMDAR model.
FUNCTION mgblock(v) {
	mgblock = 1 / (1 + exp(0.062 * -(v-sh)) * (mg / 3.57))	: from Jahr & Stevens
}
