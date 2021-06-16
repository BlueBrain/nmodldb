COMMENT

Based on the mod file written by Mark Cembrowski, 2015
Extended by a voltage-dependent outward rectification

This is an extension of the Exp2Syn class to incorporate tracking of the
specific features of different inhibitory synapses.  Specifically, this includes
whether a synapse:
	is Vgat+ (vgat)
	is Sst+ (sst)
	is Npy+ (npy)
These features are implement in order to track synapses and turn them on/off
as the simulated genotype demands.

This class also extended to have an isOn attribute, which acts as a switch
on whether the synapse is on (if = 0, conductance is always = 0; if = 1,
synapse behaves as normal).

As background, the Exp2Syn features are described as:

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
	POINT_PROCESS inhSyn
	RANGE tau1, tau2, e, i
	NONSPECIFIC_CURRENT i
	RANGE g
	RANGE vgat,sst,npy,pv,xEff,V50,slope_factor
	RANGE isOn
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	tau1=.1 (ms) <1e-9,1e9>
	tau2 = 10 (ms) <1e-9,1e9>
	e=-70	(mV)
	vgat=0
	sst=0
	npy=0
	pv=0
	xEff=-1
	isOn=0
    V50=-60.4 (mV)
    slope_factor=8.8
}

ASSIGNED {
	v (mV)
	i (nA)
	g (uS)
	factor
}

STATE {
	A (uS)
	B (uS)
}

INITIAL {
	LOCAL tp
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
	g = rect(v)*(B - A)*isOn
	i = g*(v - e)
}

DERIVATIVE state {
	A' = -A/tau1
	B' = -B/tau2
}

NET_RECEIVE(weight (uS)) {
	A = A + weight*factor
	B = B + weight*factor
}


FUNCTION rect (v(mV))( ){
	rect= 1+(0.25-1)/ ( 1. + exp (( v - V50 )/slope_factor ) ) :0.685 is the value for half activation at -60.4 mV; older model v2: 1+(0.39-1)/ ( 1. + exp (( v - -56.0 )/8.2 ) )
}
