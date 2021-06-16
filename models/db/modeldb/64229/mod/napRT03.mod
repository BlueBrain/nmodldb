NEURON { SUFFIX napRT03 }
NEURON {  USEION na READ ena WRITE ina }
ASSIGNED { ina }

PARAMETER {
	erev 		= 50    (mV)
	gmax 		= 0.4  (S/cm2)
        vrest           = 0    (mV)

	mvhalf 		= 48.
	mkconst 	= -10
	exptemp 	= 37
	mq10		= 1
	mexp 		= 1

	hvhalf 		= 0
	hkconst 	= 0
	hq10		= 1
	hexp 		= 0
        ena
} : end PARAMETER

INCLUDE "boltz_cvode.inc"

FUNCTION settau(j,v) {
  if (j==0) { : m
    if (v<-40.0) { 
      settau = .025 + .14*exp((v+40.)/10.)
    } else {
      settau = .02 + .145*exp((-v-40.)/10.)
    }
  } else {
    settau = 1
  }
}

PROCEDURE iassign () { i = g*(v-ena) ina=i }
:** kdrRT03 -- Traub kdr
