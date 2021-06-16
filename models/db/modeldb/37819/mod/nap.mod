NEURON { SUFFIX nap  }
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

	hvalence 	=  0
	hgamma		=  0
	hbaserate 	=  0
	hvhalf 		=  0
	hbasetau 	=  0
	htemp 		=  0
	hq10        =  0.
	hexp 		=  0



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
:** kdr
