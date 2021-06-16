COMMENT
------------------------------------------------------------------ 
The model here is of an NMDA synapse, but without regurds to the presynaptic
behavior. The conductance description is from: NMDA-based discrimination
in cortical neuron - Bartlett Mel. Neural Computation (1992) 4, 506.

g = gmax * ( exp(-t/tau1) - exp (-t/tau2) ) / (1 + ni *[Mg+] * exp (-gama*v) )

where gmax is 0.1- 1 nS, tau1 = 40msec, tau2 = 0.33msec, ni = 0.33 1/mM
[Mg+] =1mM, gama = 0.06 1/mV

The synaptic battery is e=0 mV


Written by M. Rapp 
Dept. of Neurobiology,
Life Science Inst.
The Hebrew University, 
Jerusalem, Israel

-----------------------------------------------------------------------
ENDCOMMENT



INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS NMDA_muki
	RANGE   g, gmax, e, i, tau1, tau2, gama, ni, Mgc, onset 
	NONSPECIFIC_CURRENT i
	 
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER {

	onset = 0 (ms)		: timing of activation
  	gmax = 0.001 (umho) 	: maximum conductance 
	tau1 = 40 (ms)
	tau2 = 0.33 (ms)
	e = 0 (mV)
	gama = 0.06 (/mv)
	ni = 0.33 (/mM)
	Mgc = 1 (mM)	: Mg+ concentration
}


ASSIGNED {
	v  (mV)		: postsynaptic voltage
	i  (nA)		: current = g*(v - Erev)
	g  (umho)		: conductance
}


 
BREAKPOINT {
	LOCAL  td
	if ( t < onset ) { 
		i = 0
	}  else {
	        td = t - onset
		g = calc_nmda ( td ) 
		i = g*(v - e)
	}
}




: This function prevents the formula to crash when V is too negetive
: and gama too large. In that case exp (-gama*v) ~= infinity.

FUNCTION calc_nmda ( td ) {
	LOCAL gama_exp

	if ( gama * v < -14 ) {  		 :implys that exp (-gama * v) > 1e6
		gama_exp = 1e6
	} else {
		gama_exp = exp(-gama*v)	  	 : the original formula
	}
	calc_nmda =  gmax * ( exp(-td/tau1) - exp(-td/tau2) )/ (1 + ni * Mgc * gama_exp )
	
}















