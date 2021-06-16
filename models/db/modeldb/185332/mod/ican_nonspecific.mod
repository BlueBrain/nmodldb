TITLE Slow Ca-dependent cation current
:
:   Ca++ dependent nonspecific cation current ICAN
:   Differential equations
:   http://senselab.med.yale.edu/modeldb/ShowModel.asp?model=144089


NEURON {
	SUFFIX ican_ns
	
	USEION ca READ cai VALENCE 2

    NONSPECIFIC_CURRENT i    
    RANGE gbar1, gbar2, i, g
    RANGE minf1, mtau1, minf2, mtau2
	GLOBAL Rb1, Rb2, caix, cac1, cac2, erev
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
}


PARAMETER {
	erev = -10	(mV)
	gbar1 = 1    (mho/cm2)
	gbar2 = 0    (mho/cm2)
    caix = 2
    cac1 = 1e-10 (mM)
    cac2 = 2e-10 (mM)

    
    
    Rb1 = 1   (/ms)
    Rb2 = 1   (/ms)
    

 
}


STATE {
	m1 m2
}

INITIAL {
    rates(cai)
    m1=minf1
    m2=minf2
}

ASSIGNED {
	i	(mA/cm2)
    
    v   (mV)	
	g       (mho/cm2)

    cai (mM)
    a1   (/ms)
    b1   (/ms)
    minf1
    mtau1    (ms)

    a2   (/ms)
    b2   (/ms)
    minf2
    mtau2    (ms)

}

BREAKPOINT { 
	SOLVE states METHOD cnexp
	g = gbar1 * m1 + gbar2 * m2

	i = g * (v - erev)
	
}

DERIVATIVE states { 
	
    rates(cai)
	m1' = (minf1 - m1) / mtau1
	m2' = (minf2 - m2) / mtau2
    
}

PROCEDURE rates(cai(mM)) {  

        

        if (cai>=0) {        
            a1 = Rb1 * (cai/cac1)^caix 
        }
        else {
            a1=0
        }

        b1 = Rb1

        mtau1 = 1/(a1+b1)
	    minf1 = a1/(a1+b1)

        if (cai>=0) {        
            a2 = Rb2 * (cai/cac2) 
        }
        else {
            a2=0
        }

        b2 = Rb2

        mtau2 = 1/(a2+b2)
	    minf2 = a2/(a2+b2)
 
}

