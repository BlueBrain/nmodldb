COMMENT

EPSP Plastic Synapse
THIS MODEL is the synaptic half of EPlasSom (EPSP Plastic Soma)

ENDCOMMENT


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
   NONSPECIFIC_CURRENT i
        POINT_PROCESS EPlasSyn
   POINTER ampa, nmda
   POINTER PreAvgCa
   POINTER postB
        RANGE precell
        RANGE gAMPA, gmaxAMPA, gNMDA, gmaxNMDA
        GLOBAL Erev_1, Erev_2
   GLOBAL AMPANMDARATIO, AMPANMDARATIO

   :SCALE
   POINTER ScaleFactor, Induction, lastprespike, AMPAMAX
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
   Erev_1 = 0           (mV)            : reversal potential
        Erev_2 = 50             (mV)            : reversal potential
   AMPANMDARATIO = 0.1

   :SCALING
   scale = 0               : USED to turn on/off scaling for some synapses during MULTI()

   :STDP
   tauLTP= 20              : time it takes to fall from 1->0.37
   tauLTD= 20              : time it takes to fall from 1->0.37
   gainLTP = 0.05           : LTP constant
   gainLTD = 0.1
   stdp = 0                : USED to turn on/off STDP
   plast = 0                : USED to accumulate plasticity during single run
   terror
   lastanyspike

}


ASSIGNED {

        dt               (ms)
        v                (mV)              : postsynaptic voltage
        i                (nA)              : current = g*(v - Erev)
        C                (mM)              : transmitter concentration
        gAMPA    (umho)    : conductance
        gNMDA    (umho) : conductance
        gmaxAMPA (umho) : maximum conductance
        gmaxNMDA (umho) : max conductance

   ampa                :PRESYNAPTIC POINTER
   nmda                :PRESYNAPTIC POINTER
   PreAvgCa            :PRESYNAPTIC POINTER
   postB               :POSTSYNAPTIC POINTER

   : SCALE
   lastprespike         :PRESYNAPTIC POINTER
   ScaleFactor          :POSTSYNAPTIC POINTER
   Induction            :POSTSYNAPTIC POINTER, used to change W at end of trial
   AMPAMAX

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
   i = gAMPA*(v-Erev_1) + gNMDA*(v-Erev_2)
}

PROCEDURE update() {

   gAMPA = gmaxAMPA * ampa
   gNMDA = gmaxNMDA * nmda * postB

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
      gmaxAMPA = gmaxAMPA*( (ScaleFactor)*scale*PreAvgCa+1 )
   }

   :STDP
   gmaxAMPA=gmaxAMPA+gmaxAMPA*plast

   :LIMITS
   if (scale > 0) {
		if (gmaxAMPA>AMPAMAX) {
			gmaxAMPA = AMPAMAX
		}
		if (gmaxAMPA<0) {
			gmaxAMPA=0
		}
	}

   :ADJUST NMDA COMPONENT
   gmaxNMDA = gmaxAMPA*AMPANMDARATIO

   :VERBATIM
    :   printf("t=%f plast=%f  W=%8.6f/%8.6f\n",t,scale,plast,gmaxNMDA,gmaxAMPA);
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
      :   printf("t=%5.2f(%5.2f/%5.2f) scale= %3.1f STDP=%3.1f  plast=%8.6f W=%8.6f/%8.6f(%5.1f)\n",t,lastprespike,lastpostspike,scale,stdp,plast,gmaxNMDA,gmaxAMPA,precell);
      :ENDVERBATIM
   }   


   VERBATIM
      return 0;
   ENDVERBATIM

}

FUNCTION STDPFunc(ISI) {
        TABLE
        DEPEND gainLTP, gainLTD
        FROM -100 TO 100 WITH 2010
        ISI=ISI-3
        if (ISI<=0.1+terror && ISI>-100) {
           STDPFunc = -exp(ISI/tauLTD)*gainLTD
        } else if (ISI>0.1+terror && ISI<100) {
           STDPFunc = exp(-ISI/tauLTP)*gainLTP
        } else {
           STDPFunc = 0
        }
}
