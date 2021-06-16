TITLE K-D channel with alpha and beta rates from RBD

UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (S) = (siemens)
        (mS) = (millisiemens)
}

NEURON {
    THREADSAFE
    : note - every variable accessible in NEURON will be having the suffix _KD2

        SUFFIX KD2
        USEION k READ ek WRITE ik
        RANGE gmax, g
}

PARAMETER {
: all values can be adjusted in hoc files
    gmax= 0.001 (S/cm2)
    np = 3		(1)

	kan=10.5	(mV)	:sets steepness of alpha
	Aan=3.0		(/ms)
	dan=-17		(mV)	:sets inflection point of beta
	kbn=30		(mV)	:sets steepness of beta
	Abn=0.12	(/ms)
	dbn=-50		(mV)
	
	kal=35		(mV)	:sets steepness of alpha
	Aal=0.00025	(/ms) :1/Aah changes tau at hyperpolarized potentials
	dal=-85		(mV)	:sets inflection point of alpha
	kbl=4.8		(mV)	:sets steepness of beta
	Abl=0.0033	(/ms) :1/Abh sets tau at depolarized potentials 
	dbl=-51		(mV)	:sets inflection point of beta
}

ASSIGNED { 
    v (mV)
    ek (mV)
    
    ik (mA/cm2)
    an (/ms)
    bn (/ms)
    al (/ms)
    bl (/ms)
	g (S/cm2)
}

STATE {
	n
	l
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    g  = gmax*n^np*l
    ik  = g*(v-ek)
}

INITIAL {
    settables(v)
    n = an/(an+bn)
    l = al/(al+bl)
}

DERIVATIVE states {  
    settables(v)      
    n' = ((an*(1-n)) - (bn*n))
    l' = ((al*(1-l)) - (bl*l))
}


PROCEDURE settables(v (mV)) {
    TABLE an, bn, al, bl DEPEND dan, dbn, dal, dbl, Aan, Abn, Aal, Abl
          FROM -100 TO 50 WITH 750

    if (v == dan) {
    	an = Aan
	} else {
		an = -Aan/kan*(dan-v) / (1 - exp((dan-v)/kan))	: forward rate of activation (/ms)
	}
	bn = Abn*exp( (dbn-v)/kbn)	: backward activation rate (/ms)
	
	al = Aal*exp( (dal-v)/kal)	: forward inactivation rate (/ms)
	bl = Abl/(exp( (dbl-v)/ kbl) + 1)	: backward inactivation rate (/ms)
}


