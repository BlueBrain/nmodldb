NEURON { SUFFIX kdrp }
NEURON { USEION k WRITE ik }
ASSIGNED { ik }
PARAMETER {
	erev 		= -90    (mV)
	gmax 		= 0.08   (umho)
        vrest           = 0

	mvalence 	= 3.
	mgamma 		= 0.7
	mbaserate 	= 10
	mvhalf 		= -30.
	mbasetau 	= 50
	mtemp 		= 24.
        mq10            = 3
	mexp 		= 2

	hvalence 	= 0
	hgamma		= 0
	hbaserate 	= 0
	hvhalf 		= 0
	hbasetau 	= 0
	htemp 		= 0
        hq10            = 3
	hexp 		= 0

	cao                (mM)
	cai                (mM)

	celsius		= 37	(degC)
	dt 				   (ms)
	v 			       (mV)

	vmax 		= 50   (mV)
	vmin 		= -100 (mV)
}
INCLUDE "bg_cvode.inc"
PROCEDURE iassign() { i = g * (v - erev) ik=i }



:** na
