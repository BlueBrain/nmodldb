COMMENT

Merged  INPUT.mod and SpkGen.mod
INPUT.mod is used as the input, and drives the synapses to the first layer Ex cells
SpkGen.mod is used to drive the input with patterns
Dean Buonomano 4/1/01

modified from pregen.mod - Mainen and Destexhe
modified to incorporate arbitrary sequences of pulses with
assigned inter-burst intervals
spkgenmode = 0 is the same as the previous version
spkgenmode = 1 permits the user to assign inter-burst intervals
in an aperiodic fashion. In this mode the on_times variable
must be initialized to +9e9 in the .oc or .hoc file.
Potential bug was fixed by adding dt/2
Dean Buonomano 5/20/99



ENDCOMMENT

UNITS {
        (mV) = (millivolt)
        (mA) = (milliamp)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    NONSPECIFIC_CURRENT i
    SUFFIX INPUT

    GLOBAL  spikedur
        RANGE   Spike, lastspk, spk

    GLOBAL seed
    GLOBAL spkgenmode
        RANGE fast_invl, slow_invl, burst_len, start, end
        RANGE noise
        RANGE on_times

}


PARAMETER {

    v
    i

    spikedur = 1.0  (ms)
        refact   = 1.5  (ms)

        spkgenmode      = 0             : 0=normal periodic mode, 1=aperiodic assigned mode
        fast_invl       = 1             : time between spikes in a burst (msec)
        slow_invl       = 100           : burst period (msec)
        burst_len       = 1             : burst length (# spikes)
        start           = 50            : location for first burst
        end             = 200           : time to stop bursting
        noise           = 0             : amount of randomness (0.0 - 1.0)
        seed            = 53892         : seed for random num generator

}

ASSIGNED {

        Spike
        spk
        lastspike

        pulsecntr
        scntr
        lcntr
        bcntr
        burst
        dt
        on_times[10]            : maximum of a sequence of 10 inter-burst intervals
                                                : should be initialized to 9e9 in .oc files

}

INITIAL {
        Spike = 0
        lastspike = -9e4

        pulsecntr = -1
        burst   = 0
        scntr   = fast_invl
        bcntr   = burst_len
        if (noise != 0) {
          set_seed(seed)
          scntr = noise * fpoisrand(fast_invl) + (1 - noise) * fast_invl
          bcntr = noise * fpoisrand(burst_len) + (1 - noise) * burst_len
        }
:       lcntr   = start - fast_invl
        lcntr   = start
}

BREAKPOINT {
        SOLVE generate
        i = Spike*(v+20)     :INSTANT STEP TO -20 mV
}



PROCEDURE generate() {
  if (t < end) {
    : leave spk up for 1 msec spike duration, Spike is not ON for ms
    if (t > lastspike+1) { spk = 0 }
    Spike = 0
    scntr = scntr - dt
    lcntr = lcntr - dt
    if (burst) {
      :  in a burst
      if (scntr <= dt+dt/2) {
                        : a spike
                        if (bcntr <= 1) {
                                : last spike in burst ?
                                burst = 0
                                if (spkgenmode==0) {
                                if (noise==0) {
                                lcntr = slow_invl
                            } else {
                            lcntr = noise * fpoisrand(slow_invl) + (1 - noise) * slow_invl
                                }
                                } else if  (spkgenmode==1) {
                                lcntr = on_times[pulsecntr] + dt
                                }
                }

                        Spike = 1
                        spk = 1
         lastspike = t

                        bcntr = bcntr - 1
                        if (noise==0) {
                                scntr = fast_invl + dt
                        } else {
                                scntr = noise * fpoisrand(fast_invl) + (1 - noise) * fast_invl
                        }
        }                                                                                       : if (scntr <= dt+dt/2)
        } else {                                                                        : end if(burst)
      :  between burst
      if (lcntr <= dt+dt/2) {                           : there is a round off problem here
                        :VERBATIM
                        :printf("FFFFFFFFFFFt = %f  pulsecntr = %f sctr = %f   x = %f   lctr=%f   spk=%f   lspk=%f\n",t,pulsecntr,scntr,x,lcntr,spk,lastspk);
                        :ENDVERBATIM
                if (noise==0) {
                        bcntr = burst_len
                } else {
                        bcntr = noise * fpoisrand(burst_len) + (1 - noise) * burst_len
                }
                burst = 1
        pulsecntr = pulsecntr + 1
      }
    }                                                   : end ~(burst)
    VERBATIM
    return 0;
    ENDVERBATIM
  }
}

FUNCTION fgauss(x,mean,std_dev) {
        fgauss = gauss(x,mean,std_dev)
}

FUNCTION fpoisrand(mean) {
  if (mean > 700) {
    fpoisrand = 4. * poisrand(mean/4.)
  } else {
    fpoisrand = poisrand(mean)
  }
}

COMMENT
Presynaptic spike generator
---------------------------

This mechanism has been written to be able to use synapses in a single
neuron receiving various types of presynaptic trains.  This is a "fake"
presynaptic compartment containing a fast spike generator.  The trains
of spikes can be either periodic or noisy (Poisson-distributed), and
either tonic or bursting.

Parameters;

   spikegenmod: 0=periodic;1=aperiodic
   noise:       between 0 (no noise-periodic) and 1 (fully noisy)
   fast_invl:   fast interval, mean time between spikes (ms)
   slow_invl:   slow interval, mean burst silent period (ms), 0=tonic train
   burst_len:   mean burst length (nb. spikes)
   on_times:    used for spkgenmode=1, assigns inter-burst intervals

Written by Z. Mainen, modified by A. Destexhe, modified by D. Buonomano

ENDCOMMENT

