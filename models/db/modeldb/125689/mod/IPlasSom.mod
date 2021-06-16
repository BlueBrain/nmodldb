COMMENT

THIS IS THE MORE EFFICIENT VERSION of EPSPPlas.
Since short-term plasticity, and the synaptic dynamics is the same
for all the presynaptic terminals of a neuron, we can calculate these
parameters presynaptically, and then pass them to each synapse model.

For computational efficiency it is obviously better to process the release
parameters in the presynaptic mechanism, rather than do the same calculations
in each synaptic mechanism.  The problem is the synaptic delay.  
To deal with this problem I have created a history (histR, histG) of the 
synaptic conductances, that are accessed by the synaptic mechanisms. The vector 
index accessed corresponds to the delay.
In the vector the first index (0) is always the current time step.
Thus if dt = 0.1 to implement a delay of 1ms you should setpointer:
setpointer syn.R_1, IN[0].soma.histR_ExIAF[10]
for a zero ms delay
setpointer syn.R_1, IN[0].soma.histR_ExIAF[0]



ENDCOMMENT


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON { 
:    NONSPECIFIC_CURRENT i
   SUFFIX IPlasSom
	RANGE C, lastrelease, lastspike, releaseat, Delay
	GLOBAL Cdur, Deadtime		
	GLOBAL Alpha, Beta 
	RANGE  gaba, R0, R1, Rinf, Rtau 
   :SHORT-TERM PLASTICITY
   GLOBAL U, trec, tfac
   RANGE  R, u, RG 

} 
 
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

STATE {
	nmda				: fraction of open NMDA channels
	R_2
}

PARAMETER {

	Cdur	= 1	  (ms)		: transmitter duration (rising phase) 
	Deadtime = 1  (ms)		: mimimum time between release events
	Delay = 0.6

	Alpha = 0.5		(/ms mM)	: 0.5 forward (binding) rate 
	Beta = 0.25		(/ms)		: 0.25 backward (unbinding) rate 
 
	:::%%% Values based on Basket Cell/type F2 from Gupta, Wang, Markram (2000)
	tfac = 20     (ms)		: this value should be close to 0 for no facilitation 
	trec = 700		 (ms) 	: recovery from depression time constant 
	U = 0.25 						: percent of transmitter released on first pulse
	
}


ASSIGNED { 
 
	dt		 (ms)
	v		 (mV)		: postsynaptic voltage
	i 		 (nA)		: current = g*(v - Erev)
	C		 (mM)		: transmitter concentration 
 
	gaba
	R0				: open channels at start of release
	R1				: open channels at end of release
	Rinf			: steady state channels open
	Rtau	(ms)		: time constant of channel binding
		 
	lastrelease	(ms)		: time of last spike 
	lastspike	(ms) 
	releaseat

	R				: Releasable pool 
	u				: for running value of U 
	
}

INITIAL {
	C = 0
	R1 = 0

	gaba = 0
	lastrelease = -9e4 
	lastspike   = -9e4 
	releaseat   = -9e4 
 
	R = 1 
	u = U 
 
}

BREAKPOINT {
    SOLVE release
:    SOLVE update
}

PROCEDURE release() { LOCAL q
    :will crash if user hasn't set pre with the connect statement
	:FIND OUT THERE WAS A SPIKE 
	q = (t - lastspike)		: time since last release ended 
						 
	if (q > Deadtime) {		: ready for another release? 
		if (v > 0) {		: spike occured? 
			lastspike = t 
			releaseat = t + Delay 
		} 
	} 
 	: CALCULATE RELEASE PARAMETERS 
	q = (t - lastrelease -Cdur)			: time since last spike with delay 
	if (q > Deadtime) {          : start release 
		if (t > releaseat - dt/2 && t < releaseat + dt/2) { 
			lastrelease = t 
     		u = u*(exptable(-q/tfac)) + U*(1-u*exptable(-q/tfac)) 
			R = R*(1 - u)*exptable(-q/trec) + 1 - exptable(-q/trec)		 
			C = R*u				: start new release, turn on 
			Rinf = C*Alpha / (C*Alpha + Beta) 
			Rtau = 1 / ((Alpha * C) + Beta) 
			R0 = gaba
		}	 						 
	} else if (q < 0) {			: still releasing? 
		: do nothing 
	} else if (C > 0) {			: in dead time after release, turn off 
		C = 0. 
		R1 = gaba 
	} 

	if (C > 0) {				: transmitter being released? 
	   gaba = Rinf + (R0 - Rinf) * exptable (- (t - lastrelease) / Rtau) 
	} else {					: no release occuring 
  	   gaba = R1 * exptable (- Beta * (t - (lastrelease + Cdur))) 
	} 

	VERBATIM
	return 0;
	ENDVERBATIM
} 


FUNCTION exptable(x) { 
	TABLE  FROM -10 TO 10 WITH 2000
 
	if ((x > -10) && (x < 10)) {
		exptable = exp(x)
	} else {
		exptable = 0.
	}
}

 
