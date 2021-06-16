TITLE Fluctuating conductances

COMMENT
-----------------------------------------------------------------------------

chris deister: This mod file was written by A. Destexhe (thank you), I
used and modified the file to approximate a noisy leak conductance. So
I got rid of inhibition (g_i) etc. resulting in only one conductance
not two. std_e and tau_e were controlled by me in the sim and usually
set to 0.001 and 5, respectively and is set in the hoc files.

20101021 ModelDB Administrator: Ted Carnevale identified three bugs
in the BREAKPOINT block.  One bug would have affected simulations 
run with tau_e == 0.  The second bug was "latent" i.e. it would 
only emerge if a naive fix were applied to the first bug--
in which case it would cause simulations run with tau_e == 0 to be 
wildly incorrect.  The third bug was that, under certain conditions, 
a minor range variable would not be updated at each fadvance.  
All three bugs have now been fixed.

Also note that the current and conductance units in this mechanism 
are all "absolute" i.e. not density units.  This is because the original 
implementation by Destexhe was as a POINT_PROCESS, but the authors are
using it here as a density mechanism without having converted to density
units.  In its present form the mechanism works in the sense that 
it can produce valid simulations.  However, its current and conductance 
values should be interpreted as having density units, i.e. mA/cm2 and S/cm2,
respectively.  In other words, the net current delivered by the mechanism 
to any compartment is the product of the numerical value of its i 
(interpreted as mA/cm2) and the compartment surface area (in cm2).  
Likewise, the net conductance presented by this mechanism to any 
compartment is the product of its ge (interpreted as S/cm2) and the 
compartment surface area.  See ModelDB for the author's
original version of syn.mod.

	Fluctuating conductance model for synaptic bombardment
	======================================================

THEORY

  Synaptic bombardment is represented by a stochastic model containing
  two fluctuating conductances g_e(t) and g_i(t) descibed by:

     Isyn = g_e(t) * [V - E_e] + g_i(t) * [V - E_i]
     d g_e / dt = -(g_e - g_e0) / tau_e + sqrt(D_e) * Ft
     d g_i / dt = -(g_i - g_i0) / tau_i + sqrt(D_i) * Ft

  where E_e, E_i are the reversal potentials, g_e0, g_i0 are the average
  conductances, tau_e, tau_i are time constants, D_e, D_i are noise diffusion
  coefficients and Ft is a gaussian white noise of unit standard deviation.

  g_e and g_i are described by an Ornstein-Uhlenbeck (OU) stochastic process
  where tau_e and tau_i represent the "correlation" (if tau_e and tau_i are 
  zero, g_e and g_i are white noise).  The estimation of OU parameters can
  be made from the power spectrum:

     S(w) =  2 * D * tau^2 / (1 + w^2 * tau^2)

  and the diffusion coeffient D is estimated from the variance:

     D = 2 * sigma^2 / tau


NUMERICAL RESOLUTION

  The numerical scheme for integration of OU processes takes advantage 
  of the fact that these processes are gaussian, which led to an exact
  update rule independent of the time step dt (see Gillespie DT, Am J Phys 
  64: 225, 1996):

     x(t+dt) = x(t) * exp(-dt/tau) + A * N(0,1)

  where A = sqrt( D*tau/2 * (1-exp(-2*dt/tau)) ) and N(0,1) is a normal
  random number (avg=0, sigma=1)


IMPLEMENTATION

  This mechanism is implemented as a nonspecific current defined as a
  point process.


PARAMETERS

  The mechanism takes the following parameters:

     E_e = -70  (mV)		: reversal potential of excitatory conductance
     E_i = -75 (mV)		: reversal potential of inhibitory conductance

     g_e0 = 0.0121 (umho)	: average excitatory conductance
     g_i0 = 0.0573 (umho)	: average inhibitory conductance

     std_e = 0.0030 (umho)	: standard dev of excitatory conductance
     std_i = 0.0066 (umho)	: standard dev of inhibitory conductance

     tau_e = 2.728 (ms)		: time constant of excitatory conductance
     tau_i = 10.49 (ms)		: time constant of inhibitory conductance



  A. Destexhe, Laval University, 1999

-----------------------------------------------------------------------------
ENDCOMMENT



INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX Gfluct
	RANGE g_e, E_e, g_e0,g_e1
	RANGE std_e, tau_e, D_e
	RANGE new_seed
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp) 
	(mV) = (millivolt)
	(umho) = (micromho)
}

PARAMETER {
	dt		(ms)

	E_e	= 0 	(mV)	: reversal potential of excitatory conductance
	g_e0	= 0.0121 (umho)	: average excitatory conductance
	std_e	= 0.0030 (umho)	: standard dev of excitatory conductance
	tau_e	= 2.728	(ms)	: time constant of excitatory conductance
}

ASSIGNED {
	v	(mV)		: membrane voltage
	i 	(nA)		: fluctuating current
	g_e	(umho)		: total excitatory conductance
	g_e1	(umho)		: fluctuating excitatory conductance
	D_e	(umho umho /ms) : excitatory diffusion coefficient
	exp_e
	amp_e	(umho)
	xtemp	(1)
}

INITIAL {
	g_e1 = 0

	if(tau_e != 0) {
		D_e = 2 * std_e * std_e / tau_e
		exp_e = exp(-dt/tau_e)
		amp_e = std_e * sqrt( (1-exp(-2*dt/tau_e)) )
	}
}

BREAKPOINT {
  SOLVE oup
  i = g_e * (v - E_e)
}

PROCEDURE oup() {
  xtemp = normrand(0,1)
  if(tau_e==0) {
    g_e1 = std_e * xtemp
    g_e = g_e1
  } else {
    g_e1 = exp_e * g_e1 + amp_e * xtemp
    g_e = g_e0 + g_e1
  }
}

PROCEDURE new_seed(seed) {		: procedure to set the seed
	set_seed(seed)
	VERBATIM
	  printf("Setting random generator with seed = %g\n", _lseed);
	ENDVERBATIM
}
