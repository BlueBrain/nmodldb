COMMENT
AMPA channel
This is an adapted version of Exp2Syn.
Adapted by Kevin M Biddell similar to as described by wolf et al 2006
4/21/07 
verified 3/29/2012 
kevin.biddell@gmail.com

Two state kinetic scheme synapse described by rise time tauon,
and decay time constant tauoff. The normalized peak condunductance is 1.
Decay time MUST be greater than rise time.

The solution of A->G->bath with rate constants 1/tauon and 1/tauoff is
 A = a*exp(-t/tauon) and
 G = a*tau2/(tauoff-tauon)*(-exp(-t/tauon) + exp(-t/tauoff))
	where tauon < tauoff

If tauoff-tauon -> 0 then we have a alphasynapse.
and if tauon -> 0 then we have just single exponential decay.

The factor is evaluated in the
initial block such that an event of weight 1 generates a
peak conductance of 1.

Because the solution is a sum of exponentials, the
coupled equations can be solved as a pair of independent equations
by the more efficient cnexp method.

ENDCOMMENT

NEURON {
	POINT_PROCESS AMPAk
	RANGE tauon, tauoff, gAmax, gA, Erev, i
	NONSPECIFIC_CURRENT i
	GLOBAL total
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
	(pS) = (picosiemens)
}

PARAMETER {
	Erev	= 0    (mV)	: reversal potential
	gAmax	= 30  (pS)	: maximal conductance fit ~5/07 by KMB
	tauon	= 1.1  (ms)<1e-9,1e9>
	tauoff	= 5.75 (ms)<1e-9,1e9>
}

ASSIGNED {
	v (mV)
	i (nA)
	gA (uS)
	factor
	total (uS)
}

STATE {
	m (uS)
	h (uS)
}

INITIAL {
	LOCAL tp
	total = 0
	if (tauon/tauoff > .9999) {
		tauon = .9999*tauoff
	}
	m = 0
	h = 0
	tp = (tauon*tauoff)/(tauoff - tauon) * log(tauoff/tauon)
	factor = -exp(-tp/tauon) + exp(-tp/tauoff)
	factor = 1/factor
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	gA = (1e-6)*gAmax*(h-m) 	: the 1e-6 is to convert pS to microSiemens
        i = gA*(v - Erev)
}	

DERIVATIVE state {
	m' = -m/tauon
	h' = -h/tauoff
}

NET_RECEIVE(weight (uS)) {
	state_discontinuity(m, m + weight*factor)
	state_discontinuity(h, h + weight*factor)
	total = total+weight
}
