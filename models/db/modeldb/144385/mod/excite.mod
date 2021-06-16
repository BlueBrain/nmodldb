TITLE excite, excitatory passive synaptic membrane channel

COMMENT
  Identical to NEURON's pas
  Renamed to place in synapses alongside pas and inhib
  Control with vector play into g.
ENDCOMMENT

UNITS {
  (mV) = (millivolt)
  (mA) = (milliamp)
  (S) = (siemens)
  (nA) = (nanoamp)
  (uS) = (microsiemens)
}
NEURON {
  SUFFIX excite
  POINT_PROCESS excite
  NONSPECIFIC_CURRENT i
  RANGE g, e
}

PARAMETER {
  g = 0.41e-3 (uS) <0,1e9> : g value will be controlled with Vector play etc.
  e = 0 (mV) : a much used excitatory reversal potential
}

ASSIGNED {v (mV) i (nA)}
BREAKPOINT {
  i = g*(v - e)
}
