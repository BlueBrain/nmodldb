NEURON {  POINT_PROCESS AMPA }
:     GLUTAMATE SYNAPSE (AMPA-Kainate receptors)
:
:     Parameters estimated from whole cell recordings of synaptic currents on
:     Cochlear neurons (Raman & Trussel, Neuron 9: 173-186, 1992) as well as
:     from sharp electrode EPSP's recordings in thalamocortical neurons (LGN)
:     (Crunelli et al. J. Physiol. 384: 603, 1987).
PARAMETER {
  Cdur	= 1.1	(ms)		: transmitter duration (rising phase)
  Alpha	= 10	(/ms mM)	: forward (binding) rate
  Beta	= 0.5	(/ms)		: backward (unbinding) rate
  Erev	= 0	(mV)		: reversal potential
  Deadtime = 2.5	(ms)		: mimimum time between release events
  GMAX  = 1	(umho)		: maximum conductance
  DELAY = 0                     : axonal delay
}
INCLUDE "sns.inc"

:* >>>> NMDA <<<<
