: Shunt current, voltage dependent
: Bob Calin-Jageman 9/20/2002
: This model implements a voltage-dependent, non-inactivating shunt current
: as described by
: Getting, 1989 and utilized by Lieb and Frost, 1998.
: 
: 
: Created by Bob Calin-Jageman
: 
: Created 	9/20/2002
: Modified 	9/20/2002
:
: Mathcheck  -  9/31/2002
: Unitscheck -  9/31/2002 (though G not literal mS)
: 
: Explanation
: Shunt current is equal to G * m * h * (V-Erev)
:	G - weight
:	m - activation level
:	h - inactivation level - always 1
:	v - current membrane potential in mv
:	Erev - reversal potential for the current
:
: m, the activation level changes as dm/dt = (m-ss - m) / tau-m
: where m-ss is the steady-state activation determined by
: m-ss = 1/( 1 + e ^((v+b)/c) )
: 
: 
: References
: 	Getting, P.A. (1989) "Reconstruction of small neural networks" in 
: Methods in Neuronal Modeling: From Synapses to Networks (1st ed), Kock & Segev
: eds, MIT Press.

: 	Lieb JR & Frost WN (1997) "Realistic Simulation of the Aplysia Siphon
: Withdrawal Reflex Circuit: Roles of Circuit Elements in Producing Motor Output"
: p. 1249 */
: 

NEURON {
	POINT_PROCESS shunt
	NONSPECIFIC_CURRENT i
	RANGE G, erev, Bm, Cm, Tm, Bh, Ch, Th, mmax, hmax, vstart
	}

UNITS {
	(mS) = (microsiemens)
	(mV) = (millivolt)
	(nA) = (nanoamp)
}

PARAMETER  {
	G = .28 (microsiemens)
	erev = -56.9 (mV)
	Bm (1)
	Cm (1)
	Tm (1)
	Bh (1)
	Ch (1)
	Th (1)
	mmax (1)
	hmax (1)
	vstart (mv)
	}

ASSIGNED {
	i	(nA)
	v	(mV)
	}

STATE { m (1) h (1)}

BREAKPOINT {
	mmax = 1/(1+exp((v+Bm)/Cm))
	hmax = 1/(1+exp((v+Bh)/Ch))
	SOLVE state METHOD cnexp
	i = G * m * h * (v - erev)
	}	

INITIAL {
	m = 1/(1+exp((vstart+Bm)/Cm))
	h = 1/(1+exp((vstart+Bh)/Ch))
}

DERIVATIVE state {
	m' = (mmax - m)/Tm
	h' = (hmax - h)/Th
}	

