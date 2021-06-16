: $Id: pregen.mod,v 1.12 1995/05/30 01:01:14 billl Exp $
: comments at end

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON	{ 
  POINT_PROCESS gen
  GLOBAL seed
  RANGE x, spk, lastspk
  RANGE fast_invl, slow_invl, burst_len, start, end
  RANGE noise
}

PARAMETER {
	fast_invl	= 1		: time between spikes in a burst (msec)
	slow_invl	= 50		: burst period (msec)
	burst_len	= 10		: burst length (# spikes)
	start		= 50		: location for first burst
	end		= 1e10		: time to stop bursting
	noise		= 0		: amount of randomness (0.0 - 1.0)
	seed            = 53892         : seed for random num generator
}

ASSIGNED {
	scntr
	lcntr
	bcntr
	burst
	x
        spk
        lastspk
	dt
}

INCLUDE "presyn.inc"

INITIAL {
	burst   = 0
	scntr	= fast_invl
	bcntr	= burst_len
	if (noise != 0) {
	  set_seed(seed)
	  scntr	= noise * fpoisrand(fast_invl) + (1 - noise) * fast_invl
	  bcntr	= noise * fpoisrand(burst_len) + (1 - noise) * burst_len
	}
	lcntr	= start - fast_invl
	x = -90
	spk = -90
}	

BREAKPOINT {
  SOLVE generate
}

PROCEDURE generate() {	

  if (t < end) {

    : leave x up for 1 msec spike duration
    if (t > lastspk+1) { spk = -90 }
    x = -90
    scntr = scntr - dt
    lcntr = lcntr - dt	

    if (burst) {
      :  in a burst

      if (scntr <= dt) {
	: a spike
	
	if (bcntr <= 1) {	
	  : last spike in burst ?

	  burst = 0
	  if (noise==0) {
	    lcntr = slow_invl
	  } else {
	    lcntr = noise * fpoisrand(slow_invl) + (1 - noise) * slow_invl
	  }
	}

	x = 50		
        spk = 50
        lastspk = t
	newspike()

	bcntr = bcntr - 1
	if (noise==0) {
	  scntr = fast_invl
	} else {
	  scntr = noise * fpoisrand(fast_invl) + (1 - noise) * fast_invl
	}
      } 	

    } else {
      :  between bursts

      if (lcntr <= dt) {
	
	if (noise==0) {
	  bcntr = burst_len
	} else {
	  bcntr = noise * fpoisrand(burst_len) + (1 - noise) * burst_len
	}
	burst = 1
	
      }

    }	

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
   noise: 	between 0 (no noise-periodic) and 1 (fully noisy)
   fast_invl: 	fast interval, mean time between spikes (ms)
   slow_invl:	slow interval, mean burst silent period (ms), 0=tonic train
   burst_len: 	mean burst length (nb. spikes)

Written by Z. Mainen, modified by A. Destexhe, The Salk Institute

ENDCOMMENT

