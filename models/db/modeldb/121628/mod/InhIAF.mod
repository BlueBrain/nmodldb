COMMENT
	
Inhibitory Integrate and Fire Unit
Buonomano 03/13/01

ENDCOMMENT

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX InhIAF 
   NONSPECIFIC_CURRENT i
	GLOBAL spikedur, refact, tauAHP, eAHP, gAHPbar
	RANGE Thr, lastspike
	RANGE gPAS, ePAS, gAHP, AHPon, gON, gOFF, eON, eOFF
}


PARAMETER {
	:v		(mv)
	gPAS = 0.001		(mho/cm2)
	ePAS = -60			(mV)

	spikedur = 1	(ms)
	refact   = 2.0	(ms)
	Thr	   = -40	(mv)

	gONconst  = 1 	(mho/cm2)
	gOFFconst = 1  (mho/cm2)
	eON  = 40		(mV)
	eOFF = -60		(mV)

	tauAHP   = 0.1 	(/ms)		: 1/tau = actual time constant of gAHP decay
	gAHPbar = 0.00005 (mho/cm2)	: peak of AHP current
	eAHP    = -90 	(mv)
}

ASSIGNED {
	v
	i		(mA/cm2)
	lastspike
	gAHP 		(mho/cm2)
	AHPon				: turns AHP on after spike ends

	gON		(mho/cm2)
	gOFF		(mho/cm2)
}


INITIAL {
	gAHP = 0
	AHPon	  = -9e4

	gON = 0
	gOFF = 0

	lastspike = -9e4
}

BREAKPOINT {
	SOLVE update
	i = gPAS*(v-ePAS) + gAHP*(v-eAHP)+gON*(v-eON)+gOFF*(v-eOFF)
}

PROCEDURE update() { LOCAL q, dv
: TURN ON AND OFF gON and gOFF to generate ACTION POTENTIAL
   gON = 0
   gOFF = 0
   q = (t-lastspike) - spikedur 

   if (q>refact) {				: refactory period over?
		if (v>Thr) {				: threshod reached?
			gON = gONconst			: turn spike current on
			lastspike = t
		}
	}
	else if ( q < 0 ) {			: spike still on
		gON=gONconst			
	} 
	else if (v > 0) {				: turn spike off
		gOFF = gOFFconst
		gAHP = gAHP + gAHPbar
		AHPon = t
	}
	gAHP = gAHP - gAHP*tauAHP*dt

	VERBATIM
		return 0;
	ENDVERBATIM

}										:END UPDATE

