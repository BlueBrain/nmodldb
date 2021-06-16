NEURON { SUFFIX kdr }
  
NEURON { USEION k WRITE ik }         

ASSIGNED { ik }

PARAMETER {
	erev 		= -75        (mV)
	gmax 		= 0.005     (umho)
        vrest           = 0     : added in with bills change to borg

	mvalence 	= 2.8
	mgamma 		=  0.7
	mbaserate 	=  .13
	mvhalf 		=  -18
	mbasetau 	=  0.3
	mtemp 		=  37
	mq10		=  3.
	mexp 		=  3

	hvalence 	= -6
	hgamma		=  0.3
	hbaserate 	=  0.095
	hvhalf 		=  -39
	hbasetau 	=  0.25
	htemp 		=  37
	hq10        =  3.
	hexp 		=  0



	cao                	 (mM)
	cai                  (mM)

	celsius			     (degC)
	dt 				     (ms)
	v 			         (mV)

	vmax 		= 100     (mV)
	vmin 		= -100   (mV)

} : end PARAMETER

INCLUDE "bg.inc"

PROCEDURE iassign () { i = g*(v-erev) ik=i }
:* SYNAPSES
:** GABAA
:________________________________________________________________
