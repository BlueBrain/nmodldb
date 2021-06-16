COMMENT

IPSP Plastic Synapse
THIS MODEL is the synaptic half of EPlasSom (EPSP Plastic Soma)

ENDCOMMENT


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON { 
   NONSPECIFIC_CURRENT i 
	POINT_PROCESS IPlasSyn
   POINTER gaba
	RANGE precell 
	RANGE gGABA, gmaxGABA
	GLOBAL Erev
	RANGE precell 
   POINTER PreAvgCa
   
   :SCALE
   POINTER ScaleFactor, Induction, lastprespike, GABAMAX 
   RANGE scale
   
   :STDP
   POINTER lastpostspike
   GLOBAL tauLTP, tauLTD, gainLTP, gainLTD
   RANGE stdp, plast, lastanyspike
   GLOBAL terror
} 
 
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}


PARAMETER {

	precell 
	Erev = -70		(mV)		: reversal potential 
   
   :SCALING            
   scale = 0               : USED to turn on/off scaling for some synapses during MULTI()

   :STDP
   tauLTP= 15              : time it takes to fall from 1->0.37
   tauLTD= 5              : time it takes to fall from 1->0.37
   gainLTP =0.8           : LTP constant
   gainLTD = 1
   stdp = 0                : USED to turn on/off STDP
   plast = 0                : USED to accumulate plasticity during single run
   terror
   lastanyspike

}


ASSIGNED { 
 
	dt		 (ms)
	v		 (mV)		   : postsynaptic voltage
	i 		 (nA)		   : current = g*(v - Erev)
	C		 (mM)		   : transmitter concentration 

   gGABA	 (umho)	   : conductance 
	gmaxGABA (umho)	: maximum conductance 
	gaba              :POINTER

   PreAvgCa            :PRESYNAPTIC POINTER
   
	: SCALE
   lastprespike         :PRESYNAPTIC POINTER
	ScaleFactor          :POSTSYNAPTIC POINTER
   Induction            :POSTSYNAPTIC POINTER, used to change W at end of trial
   GABAMAX

 	: STDP
 	lastpostspike        :POSTSYNAPTIC POINTER

   
}

INITIAL {
   terror=dt/10 
   plast = 0
   lastanyspike=-8e4
}

BREAKPOINT {
   SOLVE update
	i = gGABA*(v-Erev)
}

PROCEDURE update() {

	gGABA = gmaxGABA * gaba 
 
   if (stdp==1) {SOLVE STDP}
	if (Induction==1) {SOLVE SCALE}
	
   VERBATIM
	   return 0;
	ENDVERBATIM

} 


:%%%ONLY IS ACCESSED AT THE END OF A TRIAL%%%
PROCEDURE SCALE() { LOCAL dummy

   :SCALE
   :if (scale>0 && lastprespike>0) {
   if (scale>0) {
      :the negative number is to make IPSP decrease when ScaleFactor is Postivie
      :the constant is to change the speed by which E and IPSPs grow
      gmaxGABA = gmaxGABA*( (ScaleFactor*(-0.25))*scale*PreAvgCa+1 )
   }
   
   :STDP
   gmaxGABA=gmaxGABA+gmaxGABA*plast
   
   :LIMITS - Tiago: MAKE THIS LIMITS WORK ONLY IF THERE IS PLASTICITY. IT CAN BE ANNOYING FOR PARAMETER SEARCH TO HAVE THESE ON... DONE!
   if (scale > 0) {
		if (gmaxGABA>GABAMAX) {
			gmaxGABA = GABAMAX
		}
   if (gmaxGABA<GABAMAX/1000) {
			gmaxGABA=GABAMAX/1000
		}
	}

   :VERBATIM
   :	printf("t=%f ScaleFactor=%f  W=%8.6f/%8.6f\n",t,ScaleFactor,gmaxGABA,PreAvgCa);
   :ENDVERBATIM

   VERBATIM
   return 0;
   ENDVERBATIM

}

PROCEDURE STDP() { LOCAL dummy

   if (lastprespike-terror>lastanyspike || lastpostspike-terror>lastanyspike) {
      dummy=lastpostspike-lastprespike
      plast = plast + STDPFunc(dummy)
      :plast = plast + STDPFunc(dummy)*ScaleFactor
      if (dummy+terror>=0) {
         lastanyspike=lastpostspike
      } else if (dummy<0) {
         lastanyspike=lastprespike
      }
      :VERBATIM
      :   printf("t=%5.2f(%5.2f/%5.2f) scale= %3.1f STDP=%3.1f  plast=%8.6f W=%8.6f(%5.1f)\n",t,lastprespike,lastpostspike,scale,stdp,plast,gmaxGABA,precell);
      :ENDVERBATIM
   }   
   
   
   VERBATIM
      return 0;
   ENDVERBATIM

}

FUNCTION STDPFunc(ISI) { 
	TABLE  
	FROM -100 TO 100 WITH 2010
	: ANTI-STDP
	ISI=ISI-3
	if (ISI<=0.0+terror && ISI>-100) {
      STDPFunc = exp(ISI/tauLTP)*gainLTP
   } else if (ISI>0.0+terror && ISI<100) {
      STDPFunc = -exp(-ISI/tauLTD)*gainLTD
   } else {
      STDPFunc = 0
   }
	: STDP
	:ISI=ISI-3
	:if (ISI>0.0+terror && ISI<100) {
   :   STDPFunc = exp(ISI/tauLTP)*gainLTP
   :} else if (ISI<=0.0+terror && ISI>-100) {
   :   STDPFunc = -exp(-ISI/tauLTD)*gainLTD
   :} else {
   :   STDPFunc = 0
   :}

   : INVERTED MEXICAN HAT
	:FROM -100 TO 100 WITH 2001 
   :ISI = fabs(ISI)
   :if (ISI>=0 && ISI<100) {
   :   STDPFunc = ( -(gainLTP+gainLTD)*exp(-(ISI/tauLTD)) )+gainLTP*( exp(-(ISI/tauLTP)) ) 
   :} else {
   :   STDPFunc = 0
   :}

} 
