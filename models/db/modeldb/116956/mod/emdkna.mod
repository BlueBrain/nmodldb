COMMENT
Na dependent K channel for the fly lobular plate VS cell. This channel has a sigmoidal dependence on 
the actual Na current present.
Based on the paper:
Haag, Theunissen and Borst (1997) "The intrinsic electrophysiological characteristics of fly lobular plate tangential cells: II Active memberane properties".
J. Comp. Neurosc. 4:349-369

Author B. Torben-Nielsen @ TENU/OIST. 2009-01-13 (with help from T. Carnevale)
ENDCOMMENT

NEURON {
	SUFFIX emdkna
	USEION k READ ek WRITE ik
	USEION na READ ina
	RANGE gk, gbar, i
	RANGE ninf, ntau : would be OK for these to be GLOBAL
	GLOBAL slope,taumax
}

UNITS { : units that are not in the units database should be declared here
  (mV) = (millivolt)
  (mA) = (milliamp)
  (uA) = (microamp)
  (S) = (siemens)
}

PARAMETER {
	ek = -20 (mV) : this value will have no effect
	gbar = 0.002 (S/cm2) 	
	: midV from the paper can be left out because it is 0
	slope = 7 (uA/cm2)
	taumax = 3 (ms)
	ntau = 3 (ms) : constant
}

ASSIGNED {
	v (mV) : must declare v
	i 	(mA/cm2)
	ik 	(mA/cm2)
	gk		(S/cm2)
	ninf
	ina (mA/cm2) : for a density mechanism, ionic currents are in (mA/cm2)
}

STATE { n }

INITIAL { 
	rates(ina)
	n = ninf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gk = gbar*n*n*n*n : n^4
	i = gk * (v - ek)
	ik = i
}

DERIVATIVE states {  
	rates(ina)
	n' = (ninf - n)/ntau
}

PROCEDURE rates(ina (mA/cm2)) {
	ninf = 1 / ( 1 + exp( -ina/((1e-3)*slope) ) )  : (1e-3) converts uA/cm2 to mA/cm2
}
