NEURON {POINT_PROCESS GABAB}
PARAMETER {
  Cdur    = 150   (ms)            : transmitter duration (rising phase)
  Alpha   = 0.01  (/ms mM)        : forward (binding) rate
  Beta    = 0.005 (/ms)           : backward (unbinding) rate
  Erev    = -95   (mV)            : reversal potential (potassium)
}
INCLUDE "netcon.inc"
:* Defaults
