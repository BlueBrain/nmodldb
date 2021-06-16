NEURON { SUFFIX kmRT03 }
NEURON {  USEION k READ ek WRITE ik }
ASSIGNED { ik }

PARAMETER {
	erev 		= 0    (mV)
	gmax 		= 0    (mho/cm^2)
        vrest           = 0

	maflag 		= 2
	malphaA 	= 0.02
	malphaB		= -5.
	malphaV0	= -20.
	mbflag 		= 1
	mbetaA 		= 0.01
	mbetaB		= -18.
	mbetaV0		= -43.
	exptemp		= 37
	mq10		= 1
	mexp 		= 1

	haflag 		= 0
	halphaA 	= 0
	halphaB		= 0
	halphaV0	= 0
	hbflag 		= 0
	hbetaA 		= 0
	hbetaB		= 0
	hbetaV0		= 0
	hq10		= 3
	hexp 		= 0
        ek
} : end PARAMETER

INCLUDE "geneval_cvode.inc"

PROCEDURE iassign () { i = g*(v-ek) ik=i }
 
:** arhRT03  -- Traub ih
