NEURON { SUFFIX nap }
NEURON { USEION na WRITE ina }
ASSIGNED { ina }
PARAMETER {
	erev 		= 45        (mV)
	gmax 		= 0.1     (umho)
        vrest           = 0

	mvalence 	= 4.3
	mgamma 		=  0.7
	mbaserate 	=  4.2
	mvhalf 		=  -40
	mbasetau 	=  0.1
	mtemp 		=  37
        mq10            =  3
	mexp 		=  3

	hvalence 	=  -6
	hgamma		=  0.5
	hbaserate 	=  0.2
	hvhalf 		=  -50
	hbasetau 	=  4
	htemp 		=  37
        hq10            =  3
	hexp 		=  1

	cao                (mM)
	cai                (mM)

	celsius		= 37	     (degC)
	dt 				     (ms)
	v 			         (mV)

	vmax 		= 50     (mV)
	vmin 		= -100   (mV)
}
INCLUDE "bg_cvode.inc"
PROCEDURE iassign() { i = g * (v - erev) ina=i }
 
:* Defaults
