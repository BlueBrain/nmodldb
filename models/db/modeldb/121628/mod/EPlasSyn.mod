COMMENT

EPSP Plastic Synapse
THIS MODEL is the synaptic half of EPlasSom (EPSP Plastic Soma)

ENDCOMMENT


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON { 
    NONSPECIFIC_CURRENT i 
	POINT_PROCESS EPlasSyn
    POINTER ampa, nmda
	RANGE B, precell 
	RANGE gAMPA, gmaxAMPA, gNMDA, gmaxNMDA, Erev_1
	GLOBAL Erev_2 
    GLOBAL AMPANMDARATIO, AMPANMDARATIO, terror
    
    :SCALE
    POINTER ScaleFactor, Induction, lastprespike 
    RANGE scale    
	:LONG-TERM PLASTICITY
    POINTER lastpostspike, Gain
	RANGE LTP, Ca, stdp 	 
    GLOBAL LTPGain, LTDGain, tauLTD, tauLTP

	
} 
 
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}


PARAMETER {

	precell 
 
	Erev_1 = 0		(mV)		: reversal potential 
	Erev_2 = 50		(mV)		: reversal potential 
    AMPAMAX = 0.1               : max gmaxAMPA allowed in PLASTICITY and SCALE
    AMPANMDARATIO = 0.1
    terror
    
    :SCALING            
    scale = 1               : USED to turn on/off scaling for some synapses during MULTI()
    :STDP            
    stdp = 1
    tauLTD = 10             : decay time for LTP part of STDP function
    tauLTP = 5              : decay time for LTP part of STDP function
    LTPGain = 100
    LTDGain = 0.1
    
}


ASSIGNED { 
 
	dt		 (ms)
	v		 (mV)		: postsynaptic voltage
	i 		 (nA)		: current = g*(v - Erev)
	C		 (mM)		: transmitter concentration 
	gAMPA	 (umho)	: conductance 
	gNMDA 	 (umho)	: conductance
	gmaxAMPA (umho)	: maximum conductance 
	gmaxNMDA (umho)	: max conductance 
	B				: magnesium block
		 
	ampa                :PRESYNAPTIC POINTER
    nmda                :PRESYNAPTIC POINTER
    
	: SCALE
    lastprespike        :PRESYNAPTIC POINTER
	ScaleFactor         :POSTSYNAPTIC POINTER
    Induction           :POSTSYNAPTIC POINTER, used to change W at end of trial

   	: STDP
    lastpostspike       :POSTSYNAPTIC POINTER   
	LTP 
    Ca
	Gain                :POSTSYNAPTIC POINTER

}

INITIAL {
    terror = dt/10
	LTP = 0.0
	Ca = 0.0
}

BREAKPOINT {
    SOLVE update
	i = gAMPA*(v-Erev_1) + gNMDA*(v-Erev_2)
}

PROCEDURE update() {

    gAMPA = gmaxAMPA * ampa 
	B = mgblock(v)
	gNMDA = gmaxNMDA * nmda * B 

    if ( Induction==1 && scale==1 && lastprespike>0) {SOLVE SCALE}
    if ( stdp==1 ) {
        if ((lastprespike>t-terror && lastprespike<t+terror) ||  (lastpostspike>t-terror && lastpostspike<t+terror)) {SOLVE STDP}
    }
    VERBATIM
	    return 0;
	ENDVERBATIM
} 


PROCEDURE STDP() { LOCAL dW, dW2
					: INDUCE PLASTICITY

	dW = STDPFunc(lastpostspike-lastprespike)
	if (dW>0) {
        dW2=dW*gmaxAMPA*LTPGain*fabs(AMPAMAX-gmaxAMPA)	
	} else {
        dW2=dW*gmaxAMPA*LTDGain		
	}

    :VERBATIM
	:	printf("EPlasSyn(%3.1f) t=%6.2f, lastpre=%6.2f, lastpost=%6.2f, W=%8.6f, dW=%8.6f\n",precell,t,lastprespike,lastpostspike,gmaxAMPA,STDPFunc(lastpostspike-lastprespike));
	:ENDVERBATIM
	
	gmaxAMPA = gmaxAMPA + dW2
	
    if (gmaxAMPA>AMPAMAX) {
        gmaxAMPA = AMPAMAX
    }
    gmaxNMDA = gmaxAMPA*AMPANMDARATIO


	VERBATIM
	return 0;
	ENDVERBATIM

}


PROCEDURE SCALE() { LOCAL dummy
    
    gmaxAMPA = gmaxAMPA*ScaleFactor
    if (gmaxAMPA>AMPAMAX) {
        gmaxAMPA = AMPAMAX
    }
    gmaxNMDA = gmaxAMPA*AMPANMDARATIO
    
    :VERBATIM
    :	printf("t=%f    Induction=%f    W=%f/%f\n",t,Induction,gmaxNMDA,gmaxAMPA);
    :ENDVERBATIM

    VERBATIM
    return 0;
    ENDVERBATIM

}


 
FUNCTION mgblock(v(mV)) { 
:mgblock(-100,0,50)(w/ 0.0062 = 0.9994,0.78,0.138 
	TABLE  
	FROM -140 TO 80 WITH 1000 
	if (v>-59) { 
	   mgblock = 1/( 1+exp( (-35-v)/6 ) ) 
	} else { 
	   mgblock = 0 
	} 
} 
 

FUNCTION STDPFunc(x) { 
	TABLE  FROM -100 TO 100 WITH 2000
 
	if (x>-0.9999 && x<0.9999) {
	    STDPFunc=0
	} else if ( x>-100 && x<0)  {
	    STDPFunc = -exp(x/tauLTD)
	} else if (x>0 && x<100) {
	    STDPFunc = exp(-x/tauLTP)
	} else {
		STDPFunc = 0.
	}
}


