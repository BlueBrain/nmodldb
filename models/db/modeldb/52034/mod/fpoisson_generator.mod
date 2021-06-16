: $Id: fpoisson_generator.mod,v 1.2 2005/04/12 02:24:50 billl Exp $
TITLE Poisson_input  Poisson input generator

COMMENT

Modified for speed from poisson_generator.mod. Works if rate is not
too high (close to 1 kHz). If you need a higher rate increase no.
and/or amplitude of synapses. Paul Bush 3.95.

This input generator is designed to link up with a synapse with an interface
like the current_synapse or TrigSyn module, where a number is added to the
value of a variable to indicate that an event took place.

Usage:
object poisson_gen
poisson_gen = new fpoisson_generator(0.5)// Insert mechanism into membrane.
					// This is only necessary because of
					// neuron interface requirements.
setpointer poisson_gen.out_stim, t_synapse.stim
					// Connect to synapse.
poisson_gen.mean_rate = 0.02		// Specify rate in kHz.

ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
  POINT_PROCESS fpoisson_generator
  RANGE mean_rate, magnitude, on, off, bg_rate
  POINTER out_stim
}

ASSIGNED {
  rn
  prob1
  prob2
  dt
  out_stim
}

PARAMETER {
  mean_rate (kHz)		: Mean rate of poisson process.
  magnitude (uS)		: Magnitude of resulting EPSP.
  on = 0
  off = 1e11
  bg_rate = 0
}

INITIAL {
  prob1 = bg_rate * dt * 2^31
  prob2 = mean_rate * dt * 2^31
  :rn = (1/(2^31))
}

BREAKPOINT {
SOLVE dum
}

PROCEDURE dum() {

  out_stim = 0
  if (t >= on && t < off) {
    if (scop_random() < prob2)  {out_stim = magnitude}
  } else {
    if (scop_random() < prob1)  {out_stim = magnitude}
  }  
}
