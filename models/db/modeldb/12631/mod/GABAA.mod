NEURON {  POINT_PROCESS GABAA }
:   GABA SYNAPSE (GABA-A receptors)
:
:   Parameters estimated from whole cell recordings of GABA-A synaptic currents 
:   from dentate granule cells: Otis TS and Mody I (1992) Modulation of decay 
:   kinetics and frequency of GABA_A receptor-mediated spontaneous inhibitory
:   postsynaptic currents in hippocampal neurons.  Neurosci. 49: 13-32.
:
:   GABALOW was created to allow the use two types of GABA-A currents in the
:   same simulation.
PARAMETER {
  Cdur	= 1.0	(ms)		: transmitter duration (rising phase)
  Alpha	= 0.53	(/ms mM)	: forward (binding) rate
  Beta	= 0.18	(/ms)		: backward (unbinding) rate
  Erev	= -80	(mV)		: reversal potential
  Deadtime = 1	(ms)		: mimimum time between release events
  GMAX     = 1  (mho)		: maximum conductance
  DELAY = 0                     : axonal delay 
}
INCLUDE "sns.inc"

:* >>>> AMPA <<<<
