TITLE holt_trigsyn -- (Hopefully) faster version of TrigSyn

COMMENT
	This synapse is designed to be a very fast
synapse with an alpha function time course.  Multiple presynaptic
cells can connect to this synapse, thus obviating the need for
large numbers of synaptic mechanisms.

Usage:

First, insert this synapse in the cell at the location you are interested
in:

objectvar syn
soma {
  syn = new holt_alphasyn(x_location) // Put the synapse somewhere.
  syn.set_tau(0.1)		// Set the rise time.
}

When you want the synapse to begin, just do

syn.stim = 0.0005		// If 0.0005 uS is the maximum conductance due
				// to this input event.

ENDCOMMENT

INDEPENDENT {t FROM 0 TO 10 WITH 1 (ms)}

NEURON {
  POINT_PROCESS holt_alphasyn
  RANGE i, stim, g, erev : tau, state_t1, state_t2, scale_factor, A_coeff, B_coeff, new_state
				: Global variables.  The RANGE keyword
				: is misnamed but still necessary.
  NONSPECIFIC_CURRENT i
}

PARAMETER {
  stim = 0 (uS)			: Set this to be gmax of a synapse when you
				: want it to be triggered.
  erev = 0 (mV)			: Reversal potential for synapse.
  v (mV)
  dt (ms)
}

ASSIGNED {
  tau (ms)			: Time constant for alpha function.
  i (nA)			: Current at the present time.
  state_t1 (uS)			: Internal values used to compute g(t).
  state_t2 (uS)			: These are actually proportional to
				: g(t) at previous timesteps.
  new_state
  scale_factor			: Factor to multiply stim by to adjust state_t2
  A_coeff			: Coefficients used to compute i
  B_coeff			: (see comments below).
  g (uS)
}

INITIAL {
  IF (tau == 0.0) {
    tau = 1			: Set up a default tau.
  }
  reinit()			: Set up A and B.
}

:
: We supply a procedure to set the tau value instead of letting the user do
: it himself because various internal variables are scaled based on tau.
:
PROCEDURE set_tau(new_tau) {
  tau = new_tau			: Use the new value.
  reinit()			: Initialize internal variables.
}

FUNCTION get_tau() {		: Allow user to access the value of tau.
  get_tau = tau
}

:------------------------------------------------------------------------------
: The algorithm for computing an alpha function rapidly is based on the
: following math, demonstrated by a simple gawk program.
:
: #!/usr/local/bin/gawk -f
: #
: # Simple gawk program for producing alpha functions, to test out our
: # algorithm.
: #
: # How this works:
: #
: #  We want a simple expression relating the value of x(t) to its value at
: #  past time points.  We choose 2 past time points (since an alpha function
: #  is the product of a 2nd order differential equation) and try to construct
: #  a linear function:
: #
: #     x(t+dt) = A x(t) - B x(t-dt)
: #
: #  x(t) = t exp(-t).  Plugging this in gives:
: #
: #     (t + dt) exp(-(t+dt)) = A t exp(-t) - B (t-dt) exp(-(t-dt))
: #
: #  Solving the above equation for A and B is straightforward and it can
: #  always be done.  (Just match coefficients of t.)  The result is:
: #
: #      A = 2 exp(-dt)       B = exp(-2 dt)
: #
: #  To begin an alpha function, we pretend that it started one time step ago,
: #  i.e., we adjust x(t-dt) to be what the alpha function would have been.
: #  Since this is one time step ago, this means we simply evaluate t exp(-t)
: #  at t=-dt and subtract this result from our value for x(t-dt).
: #  
: #  This algorithm is very handy because you can make the time step as large
: #  as you like and it will still be exact.  Of course, there will be errors
: #  introduced by not sampling the function often enough, but the sampled
: #  values will be correct.
: #  
: BEGIN {
:     dt = 1;			# Time step.
:     exp_factor = exp(-dt);	# Useful constant.
:     x_t = 0;			# x(t).
:     x_tminus1 = -dt/exp_factor; # x(t-dt).  Just pretend the function started
: 				# at -dt.
:     A_coeff = 2*exp_factor;
:     B_coeff = exp_factor*exp_factor;
:     printf("# A=%g B=%g\n", A_coeff, B_coeff);
:     for (t = dt; t < 20; t += dt)
:     {
: 	x = A_coeff*x_t - B_coeff*x_tminus1;
: 	x_tminus1 = x_t;
: 	x_t = x;
: 	if (t >= 3 && t < 3+dt)	# Test starting a new superimposed spike
: 	    x_tminus1 -= dt/exp_factor;	# at time = 3.
: 	print t, x, t*exp(-t)/x;
:     }
:     exit;
: }
: -----------------------------------------------------------------------------

:
: This procedure sets up our internal variables.  It must be called every
: time the time constant or dt change.
:
PROCEDURE reinit() {
  state_t1 = 0			: Doesn't make sense to have these non-zero
  state_t2 = 0			: after changing the time constant.
  A_coeff = 2 * exp(-dt/tau)	: See above comments.
  B_coeff = exp(-2 * dt/tau)
  scale_factor = (dt/tau)*exp(dt/tau+1) : Factor to multiply state variables
				: by to set the peak of the alpha function to 1
}

:
: Main computational loop:
:
BREAKPOINT {

SOLVE dum
  g = new_state
  i = g*(v-erev)
}

PROCEDURE dum() {
    state_t2 = state_t2 - stim*scale_factor
    stim = 0
    new_state = A_coeff * state_t1 - B_coeff * state_t2
    state_t2 = state_t1		: Shift the values.
    state_t1 = new_state

VERBATIM
return 0;
ENDVERBATIM
  }
