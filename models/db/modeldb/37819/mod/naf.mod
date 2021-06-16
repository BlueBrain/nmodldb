NEURON { SUFFIX naf  }
NEURON {  USEION na READ ena WRITE ina }
ASSIGNED { ina }

PARAMETER {
	ena 		= 55       (mV)
	erev 		= 55       (mV)
	gmax 		= 0.030     (umho)

        vrest           = 0.0
	mvalence 	= 2
	mgamma 		=  0.5
	mbaserate 	=  4.5
	mvhalf 		=  -33.5
	mbasetau 	=  0.02
	mtemp 		=  36
	mq10		=  3.0
	mexp 		=  3

	hvalence 	= -6
	hgamma		=  0.3
	hbaserate 	=  0.095
	hvhalf 		=  -39
	hbasetau 	=  0.25
	htemp 		=  37
	hq10        =  3.
	hexp 		=  1



	cao                	 (mM)
	cai                  (mM)

	celsius			     (degC)
	dt 				     (ms)
	v 			         (mV)

	vmax 		=  100     (mV)
	vmin 		= -100   (mV)

} : end PARAMETER

INCLUDE "bg_cvode.inc"

PROCEDURE iassign () { i = g*(v-ena) ina=i }
:** nap
