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
	del = 1(ms)	<0,1e9>
	dur = .1(ms)	<0,1e9>
	amp (nA)
	del1 = 2(ms)	<1e-9,1e9> :time between pulses (off to on)
	n = 3		<0,1e9> : number of pulses
}
ASSIGNED { i (nA) tev (ms) cnt on a(nA)}

INITIAL {
	i = 0
	cnt = 0
	tev = del
	on = 0
	state()
}

BREAKPOINT {
	SOLVE state METHOD cvode_t
	i = a
}

PROCEDURE state() {
	if (cnt < n) {
		if (at_time(tev)) {
			if (on == 0) { : turn on
				a = amp
				on = 1
: printf("cnt=%g a=%g t=%g t-tev=%g\n", cnt+1, a, t, t-tev)
				tev = tev + dur
				at_time(tev) : notify system of the new value of tev 
			}else{ : turn off
				a = 0
				on = 0
: printf("cnt=%g a=%g t=%g t-tev=%g\n", cnt+1, a, t, t-tev)
				tev = tev + del1
                                at_time(tev) : notify system of the new value of tev 
                                cnt = cnt + 1
			}
		}
	}
}
