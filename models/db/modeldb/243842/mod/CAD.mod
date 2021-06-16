TITLE Decay of submembrane calcium concentration

: Author: Chitaranjan Mahapatra (chitaranjan@iitb.ac.in)
: Computational Neurophysiology Lab
: Indian Institute of Technology Bombay, India 

: For details refer: 
: Mahapatra C, Brain KL, Manchanda R, A biophysically constrained computational model of the action potential 
: of mouse urinary bladder smooth muscle. PLOS One (2018) 

: Inward Ca flux is maintained by T-type and L-type Ca channel


NEURON {
	SUFFIX CAD
	USEION ca READ ica,cai WRITE cai  VALENCE 2 
	RANGE depth,kd,cainf,taur, mt
}

UNITS {
	(molar) = (1/liter)			
	(mM)	= (millimolar)
	(um)	= (micron)
	(mA)	= (milliamp)
	(msM)	= (ms mM)
}

CONSTANT {
	FARADAY = 96489		(coul)		
}


PARAMETER {
	depth	= 0.41	(um)  
	taur	= 1800	(ms)		
	cainf	= 0.00015 (mM)     
	kd	= 0.0003	(mM)
	mt = -100
}

STATE {
	cai		(mM) 
}

INITIAL {
	cai = kd
}

ASSIGNED {
        ica     (mA/cm2)
	drive_channel	(mM/ms)
}
	
BREAKPOINT {
	SOLVE state METHOD derivimplicit
}

DERIVATIVE state { 

	drive_channel =  - (1000) * ica / (2 * FARADAY * depth)

	if (drive_channel <= 0.) { drive_channel = 0. }	: cannot pump below resting level

	cai' = drive_channel - mt *(cainf-cai)/taur
}

