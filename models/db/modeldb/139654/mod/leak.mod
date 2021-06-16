COMMENT
This file, leak.mod, implements the leak current
ENDCOMMENT

NEURON {
  SUFFIX leak
  NONSPECIFIC_CURRENT i
  RANGE i, Er, g
}

PARAMETER {
  g = 0 (S/cm2) < 0, 1e9 >
  Er = 0 (mV)
}

ASSIGNED {
  i (milliamp/cm2)
  v (millivolt)
}

BREAKPOINT { i = g * (v - Er) }
