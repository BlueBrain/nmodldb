TITLE inhib, inhibitory (passive) synaptic membrane channel

COMMENT
  Inhibitory synapse identical to NEURON's pas
  Renamed to place in synapses alongside pas and excite
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
  SUFFIX inhib
  POINT_PROCESS inhib
  NONSPECIFIC_CURRENT i
  RANGE g, e
}

PARAMETER {
  g = 0.41e-3 (uS) <0,1e9> : g value will be controlled with Vector play etc.
  e = -70 (mV) : an inhibitory (shunting vs hyperpol' at -80) reversal potential
}

ASSIGNED {v (mV) i (nA)}
BREAKPOINT {
  i = g*(v - e)
}
