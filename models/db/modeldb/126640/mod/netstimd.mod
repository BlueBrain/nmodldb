: $Id: netstim.mod,v 1.1.1.1 2001/01/01 20:30:37 hines Exp $
: modified in such a way that the first event will never be before start
: M.Migliore Dec.2001

NEURON	{ 
  POINT_PROCESS NetStimd
  RANGE y
  RANGE interv1, interv2, swu, swd, number, start
  RANGE noise
}

PARAMETER {
	interv1	= 10 (ms) <1e-9,1e9>: time between spikes (msec)
	interv2	= 10 (ms) <1e-9,1e9>: time between spikes (msec)
	interval	= 10
	number	= 10 <0,1e9>	: number of spikes
	start		= 50 (ms)	: start of first spike
	noise		= 0 <0,1>	: amount of randomeaness (0.0 - 1.0)
	swu		= 400
	swd		= 2000
	sw		= 2000
	flag		= 0
	tpre		= 0
}

ASSIGNED {
	y
	event (ms)
	on
	end (ms)
}

PROCEDURE seed(x) {
	set_seed(x)
}

INITIAL {
	on = 0
	y = 0
	flag=0
	tpre=0
	sw=2000
	if (noise < 0) {
		noise = 0
	}
	if (noise > 1) {
		noise = 1
	}
	if (start >= 0 && number > 0) {
		: randomize the first spike so on average it occurs at start+interval
		event = start + invl(interv1)
		net_send(event, 3)
	}
}	

PROCEDURE init_sequence(t(ms)) {
	if (number > 0) {
		on = 1
		event = t
		end = t + 1e-6 + invl(interv1)*(number-1)
	}
}

FUNCTION invl(mean (ms)) (ms) {
	if (mean <= 0.) {
		mean = .01 (ms) : I would worry if it were 0.
	}
	if (noise == 0) {
		invl = mean
	}else{
		invl = (1. - noise)*mean + noise*mean*exprand(1)
	}
}

PROCEDURE event_time() {
	if (flag==0 && t-tpre>swd) {
		flag=1 
		tpre=t
:		printf (" switch at %g %d %g \n", t, flag, swd)
	}
	if (flag==1 && t-tpre>=swu) {
		flag=0 
		tpre=t
:		printf (" switch at %g %d %g \n", t, flag, swu)
	}
	if (flag==0) {interval=interv1} else {interval=interv2}
	if (number > 0) {
		event = event + invl(interval)
	}
	if (event > end) {
		on = 0
	}
}

NET_RECEIVE (w) {
	if (flag == 0) { : external event
		if (w > 0 && on == 0) { : turn on spike sequence
			init_sequence(t)
			net_send(0, 1)
		}else if (w < 0 && on == 1) { : turn off spiking
			on = 0
		}
	}
	if (flag == 3) { : from INITIAL
		if (on == 0) {
			init_sequence(t)
			net_send(0, 1)
		}
	}
	if (flag == 1 && on == 1) {
		y = 2
		net_event(t)
		event_time()
		if (on == 1) {
			net_send(event - t, 1)
		}
		net_send(.1, 2)
	}
	if (flag == 2) {
		y = 0
	}
}

COMMENT
Presynaptic spike generator
---------------------------

This mechanism has been written to be able to use synapses in a single
neuron receiving various types of presynaptic trains.  This is a "fake"
presynaptic compartment containing a spike generator.  The trains
of spikes can be either periodic or noisy (Poisson-distributed)

Parameters;
   noise: 	between 0 (no noise-periodic) and 1 (fully noisy)
   interval: 	mean time between spikes (ms)
   number: 	mean number of spikes

Written by Z. Mainen, modified by A. Destexhe, The Salk Institute

Modified by Michael Hines for use with CVode
The intrinsic bursting parameters have been removed since
generators can stimulate other generators to create complicated bursting
patterns with independent statistics (see below)

Modified by Michael Hines to use logical event style with NET_RECEIVE
This stimulator can also be triggered by an input event.
If the stimulator is in the on=0 state and receives a positive weight
event, then the stimulator changes to the on=1 state and goes through
its entire spike sequence before changing to the on=0 state. During
that time it ignores any positive weight events. If, in the on=1 state,
the stimulator receives a negative weight event, the stimulator will
change to the off state. In the off state, it will ignore negative weight
events. A change to the on state immediately fires the first spike of
its sequence.

ENDCOMMENT

