COMMENT
This file, kleak.mod, implements the K leak current G_K(leak) in 
Quadroni and Knopfel 1994 table 2
ENDCOMMENT

NEURON {
  SUFFIX kleak
  NONSPECIFIC_CURRENT i
  RANGE i, Erev, g
}

PARAMETER {
  g = 132.8e-6 (siemens/cm2) < 0, 1e9 >
  Erev = -82 (millivolt)
}

ASSIGNED {
  i (milliamp/cm2)
  v (millivolt)
}

BREAKPOINT { i = g * (v - Erev) }
