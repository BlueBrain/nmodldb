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
   GLOBAL spikedur, refact, tauAHP, eAHP, gAHPbar, ThrConst
   RANGE Thr, lastspike
   RANGE gPAS, ePAS, gAHP, AHPon, gON, gOFF, eON, eOFF
   :FOR PLASTICITY
   RANGE SetCa, AvgCa, Ca, tauDCCa, Gain
   RANGE ScaleFactor, Induction
   GLOBAL SCALE, GainConst
   GLOBAL tstop, terror
   RANGE B                          :Mg Block for postsynaptic cell


}


PARAMETER {
        :v              (mv)
        gPAS = 0.0001			:*0.001            (mho/cm2)
        ePAS = -60                      (mV)

        spikedur = 0.6			:*1    (ms)
        refact   = 3			:*2.0  (ms)
        ThrConst = -50			:*-45
        Thr      = -50			:*-45        (mv)

        gONconst  = 1   (mho/cm2)
        gOFFconst = 1  (mho/cm2)
        eON  = 40               (mV)
        eOFF = -65				:*-60              (mV)

        tauAHP   = 0.1  (/ms)           : 0.1 1 10    1/tau = actual time constant of gAHP decay
        gAHPbar = 0.00005 (mho/cm2)     : peak of AHP current
        eAHP    = -90   (mv)

        SetCa = 1
        tauDCCa = 10                  :  in # of trials
   SCALE = 0                   : USED TO GATE SCALING PLASTICITY
   GainConst = 0.1             : LEARNING RATE
   tstop
   terror

}

ASSIGNED {
   v
   i               (mA/cm2)
   lastspike
   gAHP            (mho/cm2)
   AHPon                           : turns AHP on after spike ends

   gON             (mho/cm2)
   gOFF            (mho/cm2)

   Ca
   AvgCa
   ScaleFactor
   Induction
   Gain
   B

}


INITIAL {
   gAHP = 0
   AHPon     = -9e4

   gON = 0
   gOFF = 0
   terror = dt/10
   lastspike = -9e4

   Ca=0
   Induction = 0

}

BREAKPOINT {
        SOLVE update
        i = gPAS*(v-ePAS) + gAHP*(v-eAHP)+gON*(v-eON)+gOFF*(v-eOFF)
        B=mgblock(v)
}

PROCEDURE update() { LOCAL q, dv
: TURN ON AND OFF gON and gOFF to generate ACTION POTENTIAL
   gON = 0
   gOFF = 0
   q = (t-lastspike) - spikedur

   if (q>refact) {                              : refactory period over?
      if (v>Thr) {                            : threshod reached?
         gON = gONconst                  : turn spike current on
         lastspike = t
         Ca = Ca+1
      }
    }
    else if ( q < 0 ) {                     : spike still on
       gON=gONconst
    }
    else if (v > 0) {                               : turn spike off
       gOFF = gOFFconst
       gAHP = gAHP + gAHPbar
       AHPon = t
    }
    gAHP = gAHP - gAHP*tauAHP*dt


::: INDUCTION :::
::: HACK SO THAT INDUCTION IS RUN ON THE SECOND TO LAST TIMESTEP
    if ( t>(tstop-(2*dt)+terror) ) {
        : Induction==0(means that it has not got to INDUCTION yet
        if (Induction==0) {
            SOLVE INDUCTION
            }
        }
::: END INDUCTION :::

        VERBATIM
                return 0;
        ENDVERBATIM

}                                                                               :END UPDATE




PROCEDURE INDUCTION() {

::: INDUCTION :::
    AvgCa = AvgCa+(Ca - AvgCa)/tauDCCa
    ScaleFactor = GainConst*(SetCa-AvgCa)*SCALE
    ::: Triggers PLASTICITY IN EPSPplas
    Induction = 1
    :VERBATIM
    :     printf("INDUCTIONInhIAF             t=%f            Ca=%f           %f    %f(%f)\n",t,Ca,ScaleFactor,tstop,terror);
    :ENDVERBATIM
::: END INDUCTION :::

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
