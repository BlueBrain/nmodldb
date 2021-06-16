NEURON { SUFFIX k2RT03 }
NEURON {  USEION k READ ek WRITE ik }
ASSIGNED { ik }

PARAMETER {
  erev 		= -95    (mV)
  gmax 		= 0.4  (S/cm2)
  vrest           = 0    (mV)

  mvhalf 	= 10.
  mkconst 	= -17.
  exptemp 	= 37
  mq10		= 1
  mexp 		= 1

  hvhalf 	= 58.
  hkconst 	= 10.6
  hq10		= 1
  hexp 		= 1
  ek
} : end PARAMETER

INCLUDE "boltz_cvode.inc"

FUNCTION settau(j,v) {
  if (j==0) { : m
    settau = 4.95 + 0.5/(exp((v-81.)/25.6)+exp((-v-132.)/18.))
  } else {
    settau = 60. + 0.5/(exp((v-1.33)/200.)+exp((-v-130.)/7.1))
  }
}

PROCEDURE iassign () { i = g*(v-ek) ik=i }
 
:** kmRT03  -- Traub km
