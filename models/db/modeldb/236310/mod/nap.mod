TITLE NaP - persistent sodium current for nucleus accumbens 

COMMENT
Magistretti J, Alonso A (1999). "Biophysical properties and slow
voltage-dependent inactivation of a sustained sodium current in entorhinal
cortex layer-II principal neurons." J Gen Phys, 114: 491-509.

Traub RD, Buhl EH et al (2003). "Fast rhythmic bursting can be induced in
layer 2/3 cortical neurons by enhancing persistent na+ conductance or by
blocking BK channels." J Neurophys 89: 909-921.

Jason Moyer 2004 - jtmoyer@seas.upenn.edu
ENDCOMMENT

UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (S)  = (siemens)
}
 
NEURON {
        SUFFIX nap
        USEION na READ ena WRITE ina
        RANGE  gnabar, ina
}
 
PARAMETER {
	gnabar   =   4e-5 (S/cm2)	: 4e-5 in soma; 1.3802e-7 in dends

	mvhalf = -52.6		(mV)	: Magistretti 1999, Fig 4
	mslope = -4.6		(mV)	: Magistretti 1999, Fig 4
	mshift = 6.6594 		(mV)

	hvhalf = -48.8		(mV)	: Magistretti 1999, Fig 4
	hslope = 10.0		(mV)	: Magistretti 1999, Fig 4
	hshift = 0.52039			(mV)	

	qfact = 0.36
	
	mtau_base1 = 0.025	
	mtaux1 = 0.14
	mtau_half = -40
	mtau_slope = 10
	mtau_base2 = 0.02
	mtaux2 = 0.145
	
	
}
 
STATE { m h }
 
ASSIGNED {
	ena		(mV)
        v 		(mV)
        ina		(mA/cm2)
        gna		(S/cm2)

        minf
	hinf	

	mtau	(ms)			: Traub 2003, Table A2
   }
 
BREAKPOINT {
        SOLVE state METHOD cnexp
        gna = gnabar * m * h  
        ina = gna * ( v - ena )
:        VERBATIM
:        	printf("Ena is %g\n", ena);
:        ENDVERBATIM
}
 

 
INITIAL {
	rates(v)
	
	m = minf
	h = hinf
}

FUNCTION_TABLE tauhnap (v(mV))  (ms)		: Magistretti 1999, Fig 8A

DERIVATIVE state { 
        rates(v)
        m' = (minf - m) / mtau
        h' = (hinf - h) / (tauhnap(v)/qfact)    
}
 
PROCEDURE rates(v (mV)) {  
	TABLE minf, hinf, mtau
		FROM -200 TO 200 WITH 201

		minf = 1 / (1 + exp( (v - mvhalf - mshift) / mslope))
		hinf = 1 / (1 + exp( (v - hvhalf - hshift) / hslope))
		
		UNITSOFF
		if (v < mtau_half) {			: Traub 2003, Table A2
			mtau = mtau_base1 + mtaux1 * exp( (v - mtau_half ) / mtau_slope)
		} else {
			mtau = mtau_base2 + mtaux2 * exp( (-v + mtau_half) / mtau_slope)
		}
		UNITSON
}
 
 
