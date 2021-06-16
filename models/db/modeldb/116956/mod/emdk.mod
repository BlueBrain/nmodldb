COMMENT
K channel for the fly lobular plate VS cell. Based on the paper:
Haah, theunissen and Borst (1997) "The intrinsic electrophysiological characteristics of fly lobular plate tangential cells: II Active memberane properties".
J. Comp. Neurosc. 4:349-369

Author B. Torben-Nielsen @ TENU/OIST. 2009-01-13 (with help from T. Carnevale)
ENDCOMMENT

NEURON {
	SUFFIX emdk
	USEION k READ ek WRITE ik
	RANGE gk, gbar, i
	RANGE ninf, ntau : would be OK for these to be GLOBAL
	GLOBAL nmidv, nslope, ntaumax, nmidvdn, nslopedn, nmidvup, nslopeup
}


UNITS { : units that are not in the units database should be declared here
  (mV) = (millivolt)
  (mA) = (milliamp)
  (uA) = (microamp)
  (S) = (siemens)
}

PARAMETER {
	: set to the values described in the aforementioned paper
	ek = -20 (mV) : this value will have no effect. set in hoc code
	gbar = 0.001 (S/cm2) 	
	nmidv = 14 (mV)
	nslope = 11 (mV)
	ntaumax = 50.2 (ms)
	nmidvdn = 25 (mV)
	nslopedn = -30 (mV)
	nmidvup = 28 (mV)
	nslopeup = 32 (mV)
}

ASSIGNED {
	: either assigned by the system (e.g., v and i) or by us
	v (mV)
	i 	(mA/cm2)
	ik 	(mA/cm2)
	gk	(S/cm2)
	ninf
	ntau (ms)
}

STATE { n }

INITIAL { 
	rates(v)
	n = ninf
}

BREAKPOINT {
		SOLVE states METHOD cnexp
        gk = gbar*n*n*n*n :n^4
		i = gk * (v - ek) : for convenience, "i" is declared as range so that it can be studied as a seperate current coming from this mechanism.
		ik = i
}

DERIVATIVE states {  
		rates(v)
        n' = (ninf - n)/ntau
}

PROCEDURE rates(v (mV)) 
{
	ninf = 1/ ( 1 + exp( (nmidv-v)/nslope ) )
	ntau = ntaumax / ( exp( (nmidvdn-v)/nslopedn ) + exp( (nmidvup-v)/nslopeup ) )
	: ntau = ntau / ntaumax : EXTRA. ONLY FOR TESTING PURPOSES
}
