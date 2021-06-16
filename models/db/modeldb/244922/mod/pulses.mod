COMMENT
Since this is an electrode current, positive values of i depolarize the cell
and in the presence of the extracellular mechanism there will be a change
in vext since i is not a transmembrane current but a current injected
directly to the inside of the cell.
ENDCOMMENT

NEURON {
	POINT_PROCESS Pulses
	RANGE del, dur, amp, i, npulses, period, n
	ELECTRODE_CURRENT i
}
UNITS {
	(nA) = (nanoamp)
}

PARAMETER {
	del (ms)
	dur (ms)	<0,1e9>
	amp (nA)
	npulses
	period (ms)
}

ASSIGNED { i (nA)
		   n
		 }

INITIAL {
	i = 0
	n=0
}

BREAKPOINT {
	
	if (n==0) {
		if((t >= del) && (t<(del+dur))) { : turning on for the first time
			i = amp
		} 
		if(t >= (del+dur)){ 
			if(i==amp){:turning off for the first time
				n=n+1
			}
			i=0
		}	
	} else if (n<npulses) {
		if(t>=(del+n*period) && t<(del+n*period+dur)) {
			i=amp
		}
		if(t >= (del+n*period+dur)){ 
			if(i==amp){:turning off for the first time
				n=n+1
			}
			i=0
		}	
	}
}
