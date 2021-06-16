INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS FakeExcSyn
	POINTER pre
	RANGE lastrelease, dummy, lastspike, spike, delay, amp
	NONSPECIFIC_CURRENT i
	GLOBAL dur, Prethresh, Deadtime
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER {

	dur	= 0.1	(ms)		: current pulse duration
	amp = 2  (nA) : 2 		: pulse amplitude
	Prethresh = -20 			: voltage level nec for release
	Deadtime = 1	(ms)		: mimimum time between release events
	delay = 0 (ms)			: time between presynaptic AP and postsyn response.
							: A delay of 0.6 is used to simulate the time of AP propagation from
							: the middle of the CTX/CER axons to their synapses on TIN & RN.
}


ASSIGNED {
	i 		(nA)		: 
	pre 				: pointer to presynaptic variable
	lastrelease	(ms)		: time of last spike
	lastspike (ms)		: time of last presynaptic spike
	dummy
	spike
}

INITIAL {
	lastrelease = -9e9
	lastspike = 0
	dummy = 0
	spike = 0
}

BREAKPOINT {
	SOLVE release
	i = -dummy*amp
}

PROCEDURE release() { LOCAL q, r
	:will crash if user hasn't set pre with the connect statement 

	q = ((t - lastrelease) - dur)		: time since last release ended
	r = (t - lastspike)					: time since last spike
		
						: ready for another release?
	if (q > Deadtime && spike == 0) {
		
		if (pre > Prethresh) {		: spike occured?
			spike = 1
			lastspike = t
		}
				
	} else if (q > Deadtime && spike == 1) {
		
		if (t >= lastspike + delay) {
			dummy = 1			: start new injection
			lastrelease = t
		}
		
	} else if (q < 0) {			: still releasing?
	
		: do nothing
	
	} else if (dummy == 1) {			: in dead time after release
		dummy = 0
		spike = 0
	}
}