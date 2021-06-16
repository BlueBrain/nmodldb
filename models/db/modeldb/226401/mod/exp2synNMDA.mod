COMMENT

Combines NEURON'S Exp2Syn mechanism with 1) a previously reported
voltage-dependence of Mg2+ block (Obtained from ModelDB model #150239)
and 2) a new method to limit the total conductance to a maximum value to
produce a computationally efficient NMDAR conductance.

Explanation of the mechanism where multiple events (presynaptic
activations) cause this single NMDA synapse to saturate at most to the
supplied conductance, weight, (the blocking factor reduces the total
conductance) rather than summate:

Let's forget about the Mg2+ blocking factor for a moment and just look
at the factor g=A-B inherited from the exp2syn synapse. The g trace
rises and falls similar to an alpha function, however is formed by the
difference of two decaying exponentials A and B.  Each event resets
"g" to the rising curve before the peak (if it was already there it is
unchanged, if it was after the peak it is moved to that earlier point
that shares the same g value). Therefor at the time of an event the
value of g remains constant while the slope g' will change
from falling to rising if it is not already rising.

Tom Morse

Original comment:

Two state kinetic scheme synapse described by rise time tau1,
and decay time constant tau2. The normalized peak condunductance is 1.
Decay time MUST be greater than rise time.

The solution of A->G->bath with rate constants 1/tau1 and 1/tau2 is
 A = a*exp(-t/tau1) and
 G = a*tau2/(tau2-tau1)*(-exp(-t/tau1) + exp(-t/tau2))
	where tau1 < tau2

If tau2-tau1 -> 0 then we have a alphasynapse.
and if tau1 -> 0 then we have just single exponential decay.

The factor is evaluated in the
initial block such that an event of weight 1 generates a
peak conductance of 1.

Because the solution is a sum of exponentials, the
coupled equations can be solved as a pair of independent equations
by the more efficient cnexp method.

ENDCOMMENT

NEURON {
	POINT_PROCESS NMDA : Exp2SynNMDA
	RANGE tau1, tau2, e, i, weight, A, B
	NONSPECIFIC_CURRENT i

	RANGE g, tp, factor, cond_scale, offset
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	tau1= 0.5 (ms) : 5 (ms) .1 (ms) <1e-9, tau2>
	tau2 = 44 (ms) : 100 (ms) <tau1,1e9>
	e=0	(mV)
	mg=1    (mM)		: external magnesium concentration
        weight = 0.5e-3 (uS) : reset by events
}

ASSIGNED {
	v (mV)
	i (nA)
	g (uS)
	factor (1)
        rescale_epsps (1)
        cond_scale (1)
        tp (ms)
}

STATE {
	A (uS)
	B (uS)
}

INITIAL {
	if (tau1/tau2 > .9999) {
		tau1 = .9999*tau2
	}
	A = 0
	B = 0
	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
	factor = -exp(-tp/tau1) + exp(-tp/tau2)
	factor = 1/factor
}
BREAKPOINT {
	SOLVE state METHOD cnexp
	g = B - A
	i = g*mgblock(v)*(v - e)
}

DERIVATIVE state {
	A' = -A/tau1
	B' = -B/tau2
}

FUNCTION mgblock(v(mV)) {
	TABLE 
	DEPEND mg
	FROM -140 TO 80 WITH 1000

	: from Jahr & Stevens

	mgblock = 1 / (1 + exp(0.062 (/mV) * -v) * (mg / 3.57 (mM)))
}

: For clarity the mg voltage dependent term is ommitted from the below 
: discussion however it of course is part of the overall conductance.
: t_table returns the time on the rising curve of conductance that matches
: the conductance (so that if the conductance is falling after the peak it
: can rise to the peak again on the occurence of an event).
: To recap:  if the conductance, g, at the time of an event was falling
: the conductance is effectively moved to the same value however to
: the point on the g=A-B curve thats rising.  If the g=A-B point
: is rising when the event occurs then that point is looked up again
: in the time table, t_table, so the overall effect is that events
: occuring on the rising conductance curve have no effect (the g=A-B
: continues to rise to the peak.
: The function table values are set from hoc in set_t_table.hoc

FUNCTION_TABLE t_table(g) (ms)  

NET_RECEIVE(weight (uS)) {
        LOCAL t_tmp
        t_tmp = t_table(g/weight)
        A = weight * factor * exp( -t_tmp/tau1)
        B = weight * factor * exp( -t_tmp/tau2)
}
