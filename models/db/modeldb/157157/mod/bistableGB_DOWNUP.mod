TITLE Tristable model of AMPA receptor plasticity

COMMENT
-----------------------------------------------------------------------------
  Implementation of Graupner and Brunel's bistable model of synaptic 
  plasticity (PNAS 109(20):3991-3996, 2012.
  
  Reads internal calcium concentration (from compartment it is inserted in).
  
  Provides a state variable A that is a measure of the percentage of inserted
  AMPA receptors from fixed available pool (on scale of 0 to 1).
  
  Noise term not implemented.
  
  BPG 25-10-12
  Ausra Saudargiene 2014 04 02

-----------------------------------------------------------------------------
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX BistableGBdownup
	USEION ca READ cai
	RANGE pDOWN,pUP,pDOWN0,pUP0, C1, tau, ps, gp, gd, thp, thd, tp, td, beta, w0, w1, b, weight
	:GLOBAL tau, ps, gp, gd, thp, thd, tp, td
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(umho) = (micromho)
	(mM) = (milli/liter)
	(uM) = (micro/liter)
}

PARAMETER {
	pDOWN0	= 0 (1)	: initial synaptic efficacy
        pUP0	= 1 (1)	: initial synaptic efficacy
	tau	= 700000 (ms)	: time constant of changes
	ps	= 0.5 (1)	: transition point (unstable state)
	gp	= 1600 (1)	: potentiation rate
	gd	= 300 (1)	: depression rate
	thp	= 1.3 (uM)	: potentiation threshold re cai
	thd	= 1.0 (uM)	: depression threshold re cai
	tp	= 20000 (ms)	: time constant for LTP (eg O'Connor et al 2005)
	td	= 20000 (ms)	: time constant for LTD (eg O'Connor et al 2005)
        beta    = 0.5 (1)	: % in UP state 
	w0	= 0.5 (1)	: weight in LTD
	w1	= 1.5 (1)	: weight in LTP

}

ASSIGNED {
	cai (mM)	: calcium concentration
	C1 (uM)	        : calcium in uM
        b (1)		: ratio w1/w0
	weight (1)      : weight 
}

STATE {
	pDOWN	: synaptic efficacy (0 to 1)
        pUP	: synaptic efficacy (0 to 1)

}

INITIAL {
	pDOWN=pDOWN0
        pUP=pUP0
        b=w1/w0
        weight=(((1-wpDOWN(pDOWN))*beta+wpUP(pUP)*(1-beta))+b*(wpDOWN(pDOWN)*beta+(1-wpUP(pUP))*(1-beta)))/(beta+(1-beta)*b)

}

BREAKPOINT {
	SOLVE states METHOD euler
        
        weight=(((1-wpDOWN(pDOWN))*beta+wpUP(pUP)*(1-beta))+b*(wpDOWN(pDOWN)*beta+(1-wpUP(pUP))*(1-beta)))/(beta+(1-beta)*b)

}

DERIVATIVE states {
	
	C1 = 1000*cai

	pDOWN' = (-pDOWN*(1-pDOWN)*(ps-pDOWN)+gp*(1-pDOWN)*heavi(C1-thp)-gd*pDOWN*heavi(C1-thd))/tau

        pUP' = (-pUP*(1-pUP)*(ps-pUP)+gp*(1-pUP)*heavi(C1-thp)-gd*pUP*heavi(C1-thd))/tau

 
}

FUNCTION heavi(x (uM)) (1) {
	heavi = 0
	if (x>0) {
	  heavi = 1
	}
}

:Transition probability UP-DOWN D, argument is pUP,   D=wpUP(pUP)
FUNCTION wpUP(x) (1) {
	wpUP=0
	if (x<=0.5) {
	   wpUP = 1
	}
}

:Transition probability DOWN-UP U, argument is pDOWN, U=wpDOWN(pDOWN)
FUNCTION wpDOWN(x) (1) {
	wpDOWN=0
	if (x>=0.5) {
	   wpDOWN = 1
	}
}




