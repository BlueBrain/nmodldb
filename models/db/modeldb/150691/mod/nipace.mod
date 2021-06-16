COMMENT
Since this is an electrode current, positive values of i depolarize the cell
and in the presence of the extracellular mechanism there will be a change
in vext since i is not a transmembrane current but a current injected
directly to the inside of the cell.

n pulses of amp and dur starting at del with interval del1 between them.
ENDCOMMENT

NEURON {
	POINT_PROCESS NIClamp
	RANGE del, dur, amp, del1, n, i
	ELECTRODE_CURRENT i
}
UNITS {
	(nA) = (nanoamp)
}

PARAMETER {
	del = 100 (ms)	<0,1e9>
	dur = .1 (ms)	<0,1e9>
	amp = 2 (nA)
	del1 = 1000 (ms) <1e-9,1e9> :time between pulses (off to on)
	n = 100 	<0,1e9> : number of pulses
}
ASSIGNED { i (nA) cnt a(nA)}

INITIAL {
	i = 0
	cnt = 0
	net_send(del, 1)
}

BREAKPOINT {
	i = a
}

NET_RECEIVE(w) {
	if (flag == 1) {
		a = amp
		cnt = cnt + 1
		net_send(dur, 2)
	}else if (flag == 2) {
		a = 0
		if (cnt < n) {
			net_send(del1, 1)
		}
	}
}
