COMMENT
  IinjLT.mod
  Generates a train of photocurrent injections
  User specifies duration of each Iinj, interpulse interval (ton and toff),
  and number of pulses.
  6/30/2003 RARE LAB
           ___            ___            _ _
          /   \          /   \            |
 delay    |    \   toff  |    \           | amp
__________|     \________|     \_______  _|_
            ton
            
   num:      1              2          ...
ENDCOMMENT


NEURON {
	POINT_PROCESS IinjLTDim
	RANGE del, ton, toff, num, amp,ssI,i
	ELECTRODE_CURRENT i
}

UNITS {
	(pA) = (picoamp)
        (nA) = (nanoamp)
}

PARAMETER {
	del  = 1000 (ms)
	ton  = 8000 (ms) <0, 1e9>	: duration of ON phase
	toff = 1000 (ms) <0, 1e9>	: duration of OFF phase
	num  = 2			: how many to deliver
	amp  = 18.87043569 (pA) : absolute amplitude of the light induced current
        ssI  = 40 (pA)     : steady-state current (dark current)
        
}

ASSIGNED {
        Ncount     : counter of the number of the flashes/injections
	ival (nA)
	i (nA)
	on
	tally			: how many more to deliver
	tr (ms)   : the relative time in each flash 
        Part1
        Part2
        Part3
        ssInA (nA)
        ampnA (nA)
}

INITIAL {
	i = 0
        ssInA=0.001*ssI
        ampnA=amp*0.001
	ival = 0
	tally = num
        Ncount=0
	if (tally > 0) {
		net_send(del, 1)
		on = 0
		tally = tally - 1
	}
}

BREAKPOINT {
	LOCAL tau1,tau2,tau3,t1,t2,t3,tt,PhotoI
: printf("%g\n", t)
        tr=t-del-(ton+toff)*(Ncount-1)+dt
        if (on ==1) { 
tau1=0.140104647
tau2=0.596464897
tau3=0.590572966
t1=0.230595125
t2=0.691774291
t3=0.035251144

  
        tt=tr/1000
       	Part1=-exp(- (tt-t1 )/tau1 ) 
   	Part2= exp(- (tt-t2 )/tau2 )
   	Part3=-exp(- (tt-t3 )/tau3 )
   	PhotoI=( Part1+Part2+Part3)
       if (PhotoI>=0) {
       i=ssInA-ival*( Part1+Part2+Part3)
        } else {
       i=ssInA
	}
        
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
			ival = ampnA
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
