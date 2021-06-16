NEURON { SUFFIX kaRT03 }
NEURON {  USEION k READ ek WRITE ik }
ASSIGNED { ik }

PARAMETER {
  erev 		= -95    (mV)
  gmax 		= 0.4  (S/cm2)
  vrest           = 0    (mV)

  mvhalf 	= 60.
  mkconst 	= -8.5
  exptemp 	= 37
  mq10		= 1
  mexp 		= 4

  hvhalf 	= 78.
  hkconst 	= 6.0
  hq10		= 1
  hexp 		= 1
  ek
} : end PARAMETER

INCLUDE "boltz_cvode.inc"

FUNCTION settau(j,v) {
  if (j==0) { : m
    settau = .185 + .5/(exp((v+35.8)/19.7)+exp((-v-79.7)/12.7))
  } else {
    if (v <= -63.) {
      settau = .5/(exp((v+46.)/5.)+exp((-v-238.)/37.5))
    } else {
      settau = 9.5
    }
  }
}

PROCEDURE iassign () { i = g*(v-ek) ik=i }
:** k2RT03  -- Traub k2
