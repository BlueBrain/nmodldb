TITLE AMPA receptors with 3-state DOC kinetic scheme
COMMENT
Author: Elena Saftenku, 2003
ENDCOMMENT
NEURON {
	POINT_PROCESS AMPA_D4
      POINTER pglu
	NONSPECIFIC_CURRENT iampa
	RANGE  Erev ,Ro, RelProb
	RANGE gampa,gbarampa,ko,kd,kr,kb,kc	
	}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
	(pS) = (picosiemens)
}
PARAMETER {	
	gbarampa	= 1200   (pS)					 
	ko	= 20.35	(/ms) 				
	kc	= 4	(/ms)				
	kd	= 0.72	(/ms)					
	kr	= 0.065	(/ms)			
      kb    =0.28 (mM)  
	Erev	= 0	(mV)
      RelProb=0.46
}
ASSIGNED {
v	(mV)		
iampa (nA):current	
gampa  (pS)	: conductance 
pglu(mM): glutamate concentration
Ro : receptor occupancy	
celsius  (degC)
}
STATE {	
	C
	O
	D
}

INITIAL {
	C=1
	O=0
	D=0	
}

BREAKPOINT {
	SOLVE states METHOD derivimplicit	
	Ro=(pglu^2+ (kb^2+2*kb*pglu)*(O+D))/(pglu+kb)^2
      gampa = gbarampa * O * RelProb
	iampa = (1e-6) * gampa * (v - Erev)
}

DERIVATIVE states {
	LOCAL kon,kdn,Q10
Q10 = 2^((celsius-37(degC))/10(degC)) 
kon = ko*pglu^2/(pglu+kb)^2	
C'=Q10*(-kon*C+kc*O)
D'=Q10*(-kr*D+kd*O)
O'=Q10*(-(kc+kd)*O+kon*C+kr*D)		
}



