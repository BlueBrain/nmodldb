COMMENT
Since this is an electrode current, positive values of i depolarize the cell
and in the presence of the extracellular mechanism there will be a change
in vext since i is not a transmembrane current but a current injected
directly to the inside of the cell.
ENDCOMMENT


NEURON {
	POINT_PROCESS dbsStim
	RANGE del, dur, amp, pw, period, actPrb, adLat, active, i
	ELECTRODE_CURRENT i
}


UNITS {
	(nA) = (nanoamp)
}


PARAMETER {
	del = 0 (ms)
	dur = 1e9 (ms)	<0,1e9>
	amp = 0 (nA)
	pw = .1  (ms)
	period = 500 (ms)
	actPrb = 1
	adLat = 2
	active = 0
}


ASSIGNED { 
	i (nA) 
	offtime (ms)
}


INITIAL {
	offtime = period - pw
	i = 0
	if (del <= 0) {
		del = 0
	}
	if (active) {
		net_send(del, 1)		: initiate dbs sequence
		net_send(del+dur, 4)	: terminate dbs sequence
	}
}


BREAKPOINT {
	i = i
}


NET_RECEIVE (w) {LOCAL tmp
	if (flag == 0 && w >= 0) {
		: initiate dbs sequence from external event
		active = 1
		net_send(del, 1)
		net_send(del+dur, 4)

	} else {
	if (flag == 1 && active) {		: dbs pulse occured

		tmp = scop_random()
		if (tmp <= actPrb) {

			:may eventually want to add separate checks based on probability of propagation
			:in antidromic [net_send(adLat, 2)]
			:and orthodromic [net_event(t)] directions

			net_send(adLat, 2)		: 	turn on current after delay
			net_send(adLat+pw, 3)		:	turn off current
			net_event(t)			: signal orthodromic events to occur
		}
		net_send(adLat+offtime, 1)	: come back for next pulse

	} else {
	if (flag == 2) {
		i = amp			: up state of antidromic pulse
	} else {
	if (flag == 3) {
		i = 0				: down state of antidromic pulse
	} else {
	if (flag == 4) {
		active = 0
	}

	}
	}
	}
	}


}
