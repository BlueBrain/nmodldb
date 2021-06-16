INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON 
{
	POINT_PROCESS AMPA
	RANGE C, g, gmax, lastrelease, TRise, tau
	NONSPECIFIC_CURRENT i
	RANGE Cmax, Cdur, Alpha, Beta, Erev, Deadtime
}

UNITS 
{
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER 
{
	TRise  	= 2 (ms)
	tau    	= 2(ms)
	Cmax	= 1	(mM)		: max transmitter concentration
	Erev	= 0	(mV)		: reversal potential
	Deadtime = 1	(ms)		: mimimum time between release events
	gmax	= 0		(umho)		: maximum conductance
}


ASSIGNED 
{
	Alpha	(/ms mM)	: forward (binding) rate
	Beta	(/ms)		: backward (unbinding) rate
	Cdur	(ms)		: transmitter duration (rising phase)
	v		(mV)		: postsynaptic voltage
	i 		(nA)		: current = g*(v - Erev)
	g 		(umho)		: conductance
	C		(mM)		: transmitter concentration
	lastrelease	(ms)		: time of last spike
}

STATE
{
	R				: fraction of open channels
}

INITIAL 
{
	R = 0
	C = 0
	lastrelease = -1000
	Cdur=TRise
	Beta=1/tau
	Alpha=1/Cdur - Beta
}

BREAKPOINT 
{
	SOLVE states METHOD cnexp
	g = (gmax * R * (Alpha+Beta)) / (Alpha*(1-1/exp(1)))
	i = g*(v - Erev)
}

DERIVATIVE states
{
	evaluateC() 	: Find out value of C
	R'=Alpha * C * (1-R) - Beta * R
}

PROCEDURE evaluateC()
{
	LOCAL q
	q = ((t - lastrelease) - Cdur)		: time since last release ended
	if (q >= 0 && q <= Deadtime && C == Cmax) {	: in dead time after release
		C = 0.
	}
}

NET_RECEIVE (weight (umho)) 
{ 
	LOCAL q
	q = ((t - lastrelease) - Cdur)		: time since last release ended

: Spike has arrived, ready for another release?

	if (q > Deadtime) {
		C = Cmax			: start new release
		lastrelease = t
	} 
}

