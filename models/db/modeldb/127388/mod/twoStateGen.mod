COMMENT
Two state presynaptic event generator (burst on, burst off)
---------------------------

Generates presynaptic events at rate 1 until an external event is received.

If event is received, with probability p switches to rate 2 for a specified burst
duration, then reverts to rate 1.

Rates and duration are mean +/- std (using normal distribution for isi).

(derived from \nrn-5.8\src\nrnoc\netstim.mod)

ENDCOMMENT


NEURON	{ 
  ARTIFICIAL_CELL twoStateGen

  RANGE start
  RANGE bMean, bStd, iMean, iStd, bDurMean, bDurStd
  RANGE p
}

PARAMETER {
	start		= 0 (ms)	: start of first spike	

	bMean		= 5 (ms) 	: mean ISI during burst
	bStd		= 0
	iMean		= 100 (ms)	: mean ISI between bursts
	iStd		= 0
	bDurMean	= 50 (ms)   : mean duration of burst
	bDurStd	= 0


	p		= 1		: probability of burst onset after external event
}


ASSIGNED {
	event (ms)
	on
}


PROCEDURE seed(x) {
	set_seed(x)
}


INITIAL {
	on = 0
	if (start >= 0) {
		: randomize the first spike so on average it occurs at
		: with highest probability at start
		event = start + invl(iMean,iStd) - iMean
		: but not earlier than 0
		if (event < 0) {
			event = 0
		}
		net_send(event, 5)
	}

}


FUNCTION invl(mean (ms), stddev (ms)) {
	if (mean <= 0.) {
		mean = .01 (ms) : I would worry if it were 0.
	}

	if (stddev > 0) {
		invl = normrand(mean, stddev)
	} else {
		invl = mean
	}
}


PROCEDURE next_invl() {
	if (on == 1) {
		event = invl(bMean, bStd)
	} else {
		event = invl(iMean, iStd)
	}
	if (event < 0) {
		event = 0
	}
}


NET_RECEIVE (w) {LOCAL tmp
	if (flag == 0) { : external event
		if (w > 0 && on == 0) { : turn on burst sequence
			if (scop_random() <= p) {
				if (bDurMean > 0) {
					on = 1
					next_invl()
					net_send(event, 2)
					tmp = normrand(bDurMean, bDurStd)
					if (tmp < 0) {
						tmp = bDurStd
					}
					net_send(tmp, 3)	:shut off burst state after random duration
				}
			}

		} else {
			if (w < 0 && on == 1) { : turn off bursting
				on = 0
			}
		}

	} else {
	if (flag == 1) {
		net_event(t)
		next_invl()
		net_send(event,1)

	} else {
	if (flag == 2) {
		net_event(t)
		next_invl()
		if (on == 1) {
			net_send(event,2)		: come back for another burst spike
		}

	} else {
	if (flag == 3) {
		on = 0				: turn off burst

	} else {
	if (flag == 5) { : from INITIAL
		net_send(0,1)			: start extraburst spike sequence
	}

	}
	}
	}
	}
}
