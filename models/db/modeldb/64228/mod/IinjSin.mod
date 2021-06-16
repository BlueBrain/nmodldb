COMMENT
  IinjSin.mod
  Generates a train of sinusoidal current injections
  User specifies duration of each Iinj, interpulse interval (ton and toff),
frequency of current and total  number of sinusoidal groups.
  6/30/2003 RARE LAB
             _             _             _                _
            / \           / \           / \              / \
 delay     /   \         /   \         /   \            /   \
__________/     \       /     \       /     \__________/     \       /
            ton  \     /       \     /          toff          \     /
                  \   /         \   /                          \   /
                   \_/           \_/   ...                      \_/
                                                       
 
                        tcFreq

 number:            num1                                     num2
                                     
ENDCOMMENT


NEURON {
	POINT_PROCESS IinjSin
	RANGE del, ton, toff, num, amp,teFreq,ssI, i
	ELECTRODE_CURRENT i
}

UNITS {
	(pA) = (picoamp)
        (nA) = (nanoamp)
}

PARAMETER {
	del  = 100 (ms)
	ton  = 500 (ms) <0, 1e9>	: duration of ON phase
	toff = 1000 (ms) <0, 1e9>	: duration of OFF phase
	num  = 5			: how many to deliver
	amp  = 10 (pA)		: how big
        tcFreq = 5 (Hz)   :frequency of temporal contrast
        ssI  = 40 (pA)     : steady-state current (dark current)
}

ASSIGNED {
        Ncount     : counter of the number of the flashes/injections
	ival (nA)
	i (nA)
	on
	tally			: how many more to deliver
	tr (ms)   : the relative time in each flash 
        ssInA (nA)
}

INITIAL {
	i = 0
	ival = 0
	tally = num
        Ncount=0
        ssInA=ssI*0.001
	if (tally > 0) {
		net_send(del, 1)
		on = 0
		tally = tally - 1
	}
}

BREAKPOINT {
: printf("%g\n", t)
        tr=t-del-(ton+toff)*(Ncount-1)
        if (on ==1) { 
	i = ssInA+ival*sin(2*3.14*tcFreq*tr/1000)
        } else {
        i = ssInA+ival
        }
         
}

NET_RECEIVE (w) {
	: ignore any but self-events with flag == 1
	if (flag == 1) {
		if (on == 0) {
			: turn it on
                        Ncount=Ncount+1
			ival = amp*0.001
			on = 1
			: prepare to turn it off
			net_send(ton, 1)
		} else {
			: turn it off
			ival = 0
			on = 0
			if (tally > 0) {
				: prepare to turn it on again
				net_send(toff, 1)
				tally = tally - 1
			}
		}
	}
}
