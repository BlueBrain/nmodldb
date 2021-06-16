NEURON {
  POINT_PROCESS Shunt
  NONSPECIFIC_CURRENT i
  RANGE i, e, r
}

PARAMETER {
  r = 1 (gigaohm)  < 1e-9, 1e9 >
  e = 0 (millivolt)
}

ASSIGNED {
  i  (nanoamp)
  v  (millivolt)
}

BREAKPOINT { i = (0.001)*(v - e)/r }

