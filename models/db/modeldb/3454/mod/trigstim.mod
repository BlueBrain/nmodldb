COMMENT
Derived from IClamp but triggered by a NetCon event. On each event, a
pulse of current with dur and amp is generated.
Since this is an electrode current, positive values of i depolarize the cell
and in the presence of the extracellular mechanism there will be a change
in vext since i is not a transmembrane current but a current injected
directly to the inside of the cell.
ENDCOMMENT

NEURON {
	POINT_PROCESS TriggeredIClamp
	RANGE dur, amp, i
	ELECTRODE_CURRENT i
}
UNITS {
	(nA) = (nanoamp)
}

PARAMETER {
	dur (ms)	<0,1e9>
	amp (nA)
}

ASSIGNED { i (nA) ilocal (nA)}

INITIAL {
	i = 0
	ilocal = 0
}

BREAKPOINT {
		i = ilocal
}

NET_RECEIVE(w) {
	if (flag == 0) {
		ilocal = ilocal + amp
		net_send(dur, 1)
	}else{
		ilocal = ilocal - amp
	}
}

