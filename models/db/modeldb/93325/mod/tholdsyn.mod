: A synapse with up to 3 different conductance, each activated by a netcon event
:
: This model implements a three-conductance synapse
: This collapses into one model all three conductances
: used to describe synapses by 
: Getting, 1989 and utilized by Lieb and Frost, 1998.
: 
: Created by Bob Calin-Jageman
: 
: Created 	9/20/2002
: Modified 	9/20/2002
:
: Mathcheck  -  9/31/2002
: Unitscheck -  9/31/2002
:			  model is not meant to work in real units
:			  in particular, Go should not be interpreted as
:			  an actual conductance value in microSiemens
:			  the normalization factor also fails the units check
: 
: Explanation
: Creates a current equal to W * Go * (V-Erev) * A
: 	W - weight
:	G - current synaptic conductance
:	v - current membrane potantial
:	Erev - Reversal potential for the synapse
:	A - Normalization term
:
:	Synaptic conductance is determined by a kinetic scheme
:		tau-open - time constant of opening activated receptors
:		tau-close- time constant of closing open receptors
:
:	Pre-synaptic Action Potential ->  Gact -> Go -> Closed State
:	With every AP, Gact is set to 1
:	Activated receptors then move to Go(pen) with dGact/dt = -Gact/tau-open
:	Open receptors change as dGo/dt = Gact/tau-open - Go/tau-close
:
:	The normalization factor, A is calculated as
:	A = 1/(4e^(-3.15/(tau-close/tau-open)) + 1)
:
: Since most synapses had up to three componets, this single model
: accepts parameters for three independent conductances
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
	POINT_PROCESS tholdsyn
	NONSPECIFIC_CURRENT i
	RANGE i, G1_weight, G1_eRev, G1_opentc, G1_closetc, G2_weight, G2_eRev, G2_opentc, G2_closetc, G3_weight, G3_eRev, G3_opentc, G3_closetc, G1_A, G2_A, G3_A
	}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(S) = (microsiemens)
}

PARAMETER  {
	G1_weight = 0.0375 (1)
	G1_eRev = -80 (mV)
	G1_opentc = 10 (ms)
	G1_closetc = 25 (ms)
	G2_weight = 0.0030 (1)
	G2_eRev = -80 (mV)
	G2_opentc = 100 (ms)
	G2_closetc = 250 (ms)
	G3_weight = 0.0004 (1)
	G3_eRev = -80 (mV)
	G3_opentc = 750 (ms)
	G3_closetc = 2000 (ms)
	G1_A (1)
	G2_A (1)
	G3_A (1)
	}

ASSIGNED {
	i	(nA)
	v	(mV)
	}

STATE { G1_act (S) G1_open (S) G2_act (S) G2_open (S) G3_act (S) G3_open (S) }

BREAKPOINT {
	SOLVE states METHOD cnexp
	G1_A = (1)/((4)*(exp((-3.15)/(G1_closetc/G1_opentc)))+(1))
	G2_A = (1)/((4)*(exp((-3.15)/(G2_closetc/G2_opentc)))+(1))
	G3_A = (1)/((4)*(exp((-3.15)/(G3_closetc/G3_opentc)))+(1))
	i = (G1_weight * G1_open * (v - G1_eRev) * G1_A) + (G2_weight * G2_open * (v - G2_eRev) * G2_A) + (G3_weight * G3_open * (v - G3_eRev) * G3_A)
	}	

INITIAL {
	G1_act = 0 (S)
	G1_open = 0 (S)
	G2_act = 0 (S)
	G2_open = 0 (S)
	G3_act = 0 (S)
	G3_open = 0 (S)


}

DERIVATIVE states {
	G1_act'= -G1_act/G1_opentc
	G1_open' = G1_act/G1_opentc - G1_open/G1_closetc
	G2_act'= -G2_act/G2_opentc
	G2_open' = G2_act/G2_opentc - G2_open/G2_closetc
	G3_act'= -G3_act/G3_opentc
	G3_open' = G3_act/G3_opentc - G3_open/G3_closetc
}	

NET_RECEIVE(weight (microsiemens)) {
		state_discontinuity(G1_act, G1_act+1 (S))
		state_discontinuity(G2_act, G2_act+1 (S))
		state_discontinuity(G3_act, G3_act+1 (S))
}