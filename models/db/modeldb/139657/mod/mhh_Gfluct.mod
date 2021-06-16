TITLE Fluctuating conductances

COMMENT
-----------------------------------------------------------------------------

	Fluctuating conductance model for synaptic bombardment
	======================================================
	
This implementation models only excitatory synaptic input, and is based on the model described in:

Rudolph M, Destexhe A (2005) An extended analytic expression for the membrane potential distribution of conductance-based synaptic noise. Neural Comput 17:2301-15

Original code:
http://senselab.med.yale.edu/ModelDb/showmodel.asp?model=64259&file=\NCnote\Gfluct.mod
  

IMPLEMENTATION

  This mechanism is implemented as a nonspecific current defined as a
  point process.


PARAMETERS

  The mechanism takes the following parameters:

     E_e (mV)		: reversal potential of excitatory conductance

     g_e0 (umho)	: average excitatory conductance

     std_e (umho)	: standard dev of excitatory conductance

     tau_e (ms)		: time constant of excitatory conductance

  A. Destexhe, Laval University, 1999

  Trivial modifications by:
  Matthias H. Hennig, University of Edinburgh, 2011
  
-----------------------------------------------------------------------------
ENDCOMMENT



INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS mhh_Gfluct
	RANGE g_e, E_e, g_e0, g_e1
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
	g_e0	= 0.000001 (umho)	: average excitatory conductance
	std_e	= 0.0002 (umho)	: standard dev of excitatory conductance
	tau_e	= 2	(ms)	: time constant of excitatory conductance
}

ASSIGNED {
	v	(mV)		: membrane voltage
	i 	(nA)		: fluctuating current
	g_e	(umho)		: total excitatory conductance
	g_e1	(umho)		: fluctuating excitatory conductance
	D_e	(umho umho /ms) : excitatory diffusion coefficient
	exp_e
	amp_e	(umho)
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
	if(tau_e==0) {
	   g_e = std_e * normrand(0,1)
	}
	g_e = g_e0 + g_e1
	i = g_e * (v - E_e)
}


PROCEDURE oup() {		: use Scop function normrand(mean, std_dev)
   if(tau_e!=0) {
	g_e1 =  exp_e * g_e1 + amp_e * normrand(0,1)
   }
}


PROCEDURE new_seed(seed) {		: procedure to set the seed
	set_seed(seed)
	VERBATIM
	  printf("Setting random generator with seed = %g\n", _lseed);
	ENDVERBATIM
}

