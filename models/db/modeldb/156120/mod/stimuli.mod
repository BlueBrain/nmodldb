COMMENT
Since this is an electrode current, positive values of i depolarize the cell
and in the presence of the extracellular mechanism there will be a change
in vext since i is not a transmembrane current but a current injected
directly to the inside of the cell.
ENDCOMMENT

NEURON {
	POINT_PROCESS IC_altern
	RANGE del, dur, amp, i,freq
	ELECTRODE_CURRENT i
}
UNITS {
	(nA) = (nanoamp)
	 PI  = (pi) (1)
}

PARAMETER {
	del (ms)
	dur (ms)	<0,1e9>
	amp (nA)
	freq (1/s)

}
ASSIGNED { i (nA) }

INITIAL {
	i = 0

}

BREAKPOINT {
	at_time(del)
	at_time(del+dur)

	if (t >= del&&t < del + dur) {
		i = amp*sin(2*PI*freq*(t-del)*(1e-3))
	}else{
		i = 0
	}
}
