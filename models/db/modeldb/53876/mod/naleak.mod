COMMENT
This file, naleak.mod, implements the na leak current G_Na(leak) in 
Quadroni and Knopfel 1994 table 2
ENDCOMMENT

NEURON {
  SUFFIX naleak
  NONSPECIFIC_CURRENT i
  RANGE i, Erev, g
}

PARAMETER {
  g = 30.2e-6 (siemens/cm2) < 0, 1e9 >
  Erev = 50 (millivolt)
}

ASSIGNED {
  i (milliamp/cm2)
  v (millivolt)
}

BREAKPOINT { i = g * (v - Erev) }
