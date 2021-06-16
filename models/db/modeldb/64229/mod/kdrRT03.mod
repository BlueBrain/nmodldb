NEURON { SUFFIX kdrRT03 }
NEURON {  USEION k READ ek WRITE ik }
ASSIGNED { ik }

PARAMETER {
  erev 		= -95    (mV)
  gmax 		= 0.4  (S/cm2)
  vrest           = 0    (mV)

  mvhalf 	= 29.5
  mkconst 	= -10
  exptemp 	= 37
  mq10		= 1
  mexp 		= 4

  hvhalf 	= 0
  hkconst 	= 0
  hq10		= 1
  hexp 		= 0
  ek
} : end PARAMETER

INCLUDE "boltz_cvode.inc"

FUNCTION settau(j,v) {
  if (j==0) { : m
    if (v<-10.0) { 
      settau = .25 + 4.35*exp((v+10.)/10.)
    } else {
      settau = .25 + 4.35*exp((-v-10.)/10.)
    }
  } else {
    settau = 1
  }
}

PROCEDURE iassign () { i = g*(v-ek) ik=i }
 
:** kaRT03  -- Traub ka
