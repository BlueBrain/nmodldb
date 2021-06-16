: Leak Current Model, Bob Calin-Jageman 9/20/2002
: This model implements a simple leak function described by
: Getting, 1989 and utilized by Lieb and Frost, 1998.
: 
: 
: Created by Bob Calin-Jageman
: 
: Created 	9/20/2002
: Modified 	9/20/2002
: 
: Mathcheck  - Complete  9/31/2002
: Unitscheck - questions about resistance units
: 
: Explanation
: Leak current calculated as difference between current and resting membrane
: potentials divided by resistance
:	
: Resistance should be in megaOhms
: resting potential should be in millivolts
: output current will be in milliamps
: 
: References
: 	Getting, P.A. (1989) "Reconstruction of small neural networks" in 
: Methods in Neuronal Modeling: From Synapses to Networks (1st ed), Kock & Segev
: eds, MIT Press.

: 	Lieb JR & Frost WN (1997) "Realistic Simulation of the Aplysia Siphon
: Withdrawal Reflex Circuit: Roles of Circuit Elements in Producing Motor Output"
: p. 1249 */
: 
NEURON {
	POINT_PROCESS leak
	NONSPECIFIC_CURRENT i
	RANGE r, vrest
	}

UNITS {
	(mV) = (millivolt)
	(nA) = (nanoamp)
}

PARAMETER  {
	r = 15.7 (megaohm)
	vrest = -45 (mV)
	}

ASSIGNED {
	i	(nA)
	v	(mV)
	}


BREAKPOINT {
	i = (v - vrest)/r
	}	
