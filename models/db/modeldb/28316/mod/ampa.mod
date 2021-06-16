TITLE simple AMPA receptors

COMMENT
-----------------------------------------------------------------------------

	Simple model for glutamate AMPA receptors
	=========================================

  - FIRST-ORDER KINETICS, FIT TO WHOLE-CELL RECORDINGS

    Whole-cell recorded postsynaptic currents mediated by AMPA/Kainate
    receptors (Xiang et al., J. Neurophysiol. 71: 2552-2556, 1994) were used
    to estimate the parameters of the present model; the fit was performed
    using a simplex algorithm (see Destexhe et al., J. Computational Neurosci.
    1: 195-230, 1994).

  - SHORT PULSES OF TRANSMITTER (0.3 ms, 0.5 mM)

    The simplified model was obtained from a detailed synaptic model that 
    included the release of transmitter in adjacent terminals, its lateral 
    diffusion and uptake, and its binding on postsynaptic receptors (Destexhe
    and Sejnowski, 1995).  Short pulses of transmitter with first-order
    kinetics were found to be the best fast alternative to represent the more
    detailed models.

  - ANALYTIC EXPRESSION

    The first-order model can be solved analytically, leading to a very fast
    mechanism for simulating synapses, since no differential equation must be
    solved (see references below).



References

   Destexhe, A., Mainen, Z.F. and Sejnowski, T.J.  An efficient method for
   computing synaptic conductances based on a kinetic model of receptor binding
   Neural Computation 6: 10-14, 1994.  

   Destexhe, A., Mainen, Z.F. and Sejnowski, T.J. Synthesis of models for
   excitable membranes, synaptic transmission and neuromodulation using a 
   common kinetic formalism, Journal of Computational Neuroscience 1: 
   195-230, 1994.

-----------------------------------------------------------------------------
ENDCOMMENT



NEURON {
	POINT_PROCESS AMPA_S
	NONSPECIFIC_CURRENT i
	RANGE R, g, gmax, i
	GLOBAL Cdur_a, Alpha_a, Beta_a, Erev_a, Rinf_a, Rtau_a
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER {

	Cdur_a	= 1	(ms)		: transmitter duration (rising phase)
	Alpha_a	= 1.1	(/ms)	: forward (binding) rate
	Beta_a	= 0.19	(/ms)		: backward (unbinding) rate
	Erev_a	= 0	(mV)		: reversal potential
	gmax 
}


ASSIGNED {
	v		(mV)		: postsynaptic voltage
	i		(nA)		: current = g*(v - Erev)
	g 		(umho)		: conductance
	Rinf_a				: steady state channels open
	Rtau_a		(ms)		: time constant of channel binding
	synon
}

STATE {Ron Roff}

INITIAL {
	Rinf_a = Alpha_a / (Alpha_a + Beta_a)
	Rtau_a = 1 / (Alpha_a + Beta_a)
	synon = 0
}

BREAKPOINT {
	SOLVE release METHOD cnexp
	g = gmax*(Ron + Roff)*1(umho)
	i = g*(v - Erev_a)
}

DERIVATIVE release {
	Ron' = (synon*Rinf_a - Ron)/Rtau_a
	Roff' = -Beta_a*Roff
}

: following supports both saturation from single input and
: summation from multiple inputs
: if spike occurs during CDur then new off time is t + CDur
: ie. transmitter concatenates but does not summate
: Note: automatic initialization of all reference args to 0 except first

NET_RECEIVE(weight, on, nspike, r0, t0 (ms)) {
	: flag is an implicit argument of NET_RECEIVE and  normally 0
        if (flag == 0) { : a spike, so turn on if not already in a Cdur_a pulse
		nspike = nspike + 1
		if (!on) {
			r0 = r0*exp(-Beta_a*(t - t0))
			t0 = t
			on = 1
			synon = synon + weight
			state_discontinuity(Ron, Ron + r0)
			state_discontinuity(Roff, Roff - r0)
		}
		: come again in Cdur_a with flag = current value of nspike
		net_send(Cdur_a, nspike)
        }
	if (flag == nspike) { : if this associated with last spike then turn off
		r0 = weight*Rinf_a + (r0 - weight*Rinf_a)*exp(-(t - t0)/Rtau_a)
		t0 = t
		synon = synon - weight
		state_discontinuity(Ron, Ron - r0)
		state_discontinuity(Roff, Roff + r0)
		on = 0
	}
}

