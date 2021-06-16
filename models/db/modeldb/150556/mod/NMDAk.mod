COMMENT
NMDA channel
This is an adapted version of Exp2Syn.
Adapted by Kevin M Biddell similar to as described by wolf et al 2006
4/21/07
verified 3/29/2012 KMB 
kevin.biddell@gmail.com

Voltage dependence of Mg2+ block:
	Jahr & Stevens 1990. J Neurosci 10: 1830.
	Jahr & Stevens 1990. J Neurosci 10: 3178.

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
	POINT_PROCESS NMDAk
	RANGE tauon, tauoff, gNmax, gN, Erev, i,mg
	NONSPECIFIC_CURRENT i
	GLOBAL total,vmin, vmax
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
	(pS) = (picosiemens)
	(mM) = (milli/liter)
}

PARAMETER {
	Erev	= 0    (mV)	: reversal potential
	gNmax	= 60  (pS)	: maximal conductance :fit by KMB ~5/07
	tauon	= 2.23  (ms)<1e-9,1e9>	: changed by KMB 5/07/07 from 2.82 to 20.77% less (wolf et 2006)
	tauoff	= 75.68 (ms)<1e-9,1e9> 	: changed by KMB 4/23/07 from 160 to 52.7% less (wolf et al 2006)
	mg	= 1.2    (mM)	: external magnesium concentration
	vmin = -120	(mV)
	vmax = 100	(mV)
}

ASSIGNED {
	v (mV)
	i (nA)
	gN (uS)
	factor
	total (uS)
}

STATE {
	m (uS)
	h (uS)
	B		: fraction free of Mg2+ block
}

INITIAL {
	LOCAL tp
	rates(v)
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
	rates(v)
	SOLVE state METHOD cnexp
	gN = (1e-6)*gNmax*(h-m) 	: the 1e-6 is to convert pS to microSiemens
        i = gN*B*(v - Erev)
	
}

DERIVATIVE state {
	m' = -m/tauon
	h' = -h/tauoff
}

PROCEDURE rates(v(mV)) {
	TABLE B
	DEPEND mg
	FROM vmin TO vmax WITH 200

	: from Jahr & Stevens

	B = 1 / (1 + exp(0.062 (/mV) * -v) * (mg / 3.57 (mM)))
}


NET_RECEIVE(weight (uS)) {
	state_discontinuity(m, m + weight*factor)
	state_discontinuity(h, h + weight*factor)
	total = total+weight
}
