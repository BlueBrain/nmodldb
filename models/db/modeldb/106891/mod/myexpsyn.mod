: $Id: myexpsyn.mod,v 1.2 2007/12/28 02:51:31 billl Exp $

COMMENT
Functionally unchanged from original expsyn.mod which is POINT_PROCESS ExpSyn
Added are RANGE variables poid,poty,prid,prty which are
  postsynaptic_ID, presynaptic_ID, postsynaptic_cell_type, presynaptic_cell_type
  respectively
This was used to insure that JitCon hooks itself up to the correct synapse -- this turns out
  not to be necessary since the correct synapse will always be the one at the end of the list
ENDCOMMENT

NEURON {
	POINT_PROCESS expsyn
	RANGE tau, e, i, poid, poty, prid, prty
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	tau = 0.1 (ms) <1e-9,1e9>
	e = 0	(mV)
}

ASSIGNED {
	v (mV)
	i (nA)
        poid poty prid prty
}

STATE {
	g (uS)
}

INITIAL {
	g=0
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	i = g*(v - e)
}

DERIVATIVE state {
	g' = -g/tau
}

NET_RECEIVE(weight (uS)) {
	state_discontinuity(g, g + weight)
}
