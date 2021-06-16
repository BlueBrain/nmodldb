TITLE NMDA synapse

COMMENT
	simple alpha-synapse that generates a single PSP 
        *********************************************
        reference:      McCormick, Wang & Huguenard (1993)
                        Cerebral Cortex 3(5), 387-398
        *********************************************
	Assembled for MyFirstNEURON by Arthur Houweling
ENDCOMMENT
					       
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS NmdaSynapse
	USEION mg READ mgo
	RANGE onset, gmaxEPSP, e, g, i, w
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) 	= (nanoamp)
	(mV)	= (millivolt)
	(nS) 	= (nanomho)
	(mM)    = (milli/liter)
        F	= 96480 (coul)
        R       = 8.314 (volt-coul/degC)

}

PARAMETER {
	onset= 20	(ms)
	gmaxEPSP= 0	(nS)
	w= 1				: weight factor for gmaxEPSP
	e= 0.0		(mV)
	v		(mV)
	mgo		(mM)
	celsius 	(degC)
}

ASSIGNED { 
	i (nA)  
	g (nS)
	tadj
}

UNITSOFF
INITIAL {
        tadj = 3^((celsius-23.5)/10)
}    

BREAKPOINT {  LOCAL k, tt
	tt= (t-5-onset)*tadj
	k = 1.07 * exp(2*0.73*v*(.001)*F/(R*(celsius+273.16)))
	if (t>onset+5) {
	  g = w*gmaxEPSP * (1-(1/(1+k/mgo))) * 
                            exp(-tt/71) * (1-exp(-tt/13.2))/0.60
	}
	else { g = 0}
	i = g * (v-e)
}
UNITSON

