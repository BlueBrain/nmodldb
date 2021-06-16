: $Id: pregen.mod,v 1.3 2000/05/16 11:16:56 hines Exp $
: comments at end

NEURON	{ 
  POINT_PROCESS SpikeGenerator2
  RANGE y
  RANGE start, end,on
  RANGE time
}

PARAMETER {
	start		= 100 (ms)	: start of first interburst interval
	end		= 1e10 (ms)	: time to stop bursting
}

ASSIGNED {
	y
	 
	event (ms)
	 

	on

	time[200] (ms)	
	indice
}


INITIAL {	
	 
	indice=0
	on = 1
	y = -90
	event = start 
	event_time()
	while (on == 1 && event < 0) {
		event_time()
	}
	if (on == 1) {
		net_send(event, 1)
	}
}	

PROCEDURE event_time() {
	event=event+time[indice] 
	indice=indice+1
	if (event > end) {
		on = 0
	}
}

NET_RECEIVE (w) {
	if (flag == 1 && on == 1) {
		y = 20
		net_event(t)
		event_time()
		net_send(event - t, 1)
		net_send(.1, 2)
	}
	if (flag == 2) {
		y = -90
	}
}

COMMENT
Presynaptic spike generator
---------------------------

This mechanism has been written to be able to use synapses in a single
neuron receiving various types of presynaptic trains.  This is a "fake"
presynaptic compartment containing a fast spike generator.  The trains
of spikes can be either periodic or noisy (Poisson-distributed), and 
either tonic or bursting.

Parameters;
   noise: 	between 0 (no noise-periodic) and 1 (fully noisy)
   fast_invl: 	fast interval, mean time between spikes (ms)
   slow_invl:	slow interval, mean burst silent period (ms), 0=tonic train
   burst_len: 	mean burst length (nb. spikes)

Written by Z. Mainen, modified by A. Destexhe, The Salk Institute

Modified by Michael Hines for use with CVode

Modified by Michael Hines to use logical event style with NET_RECEIVE
ENDCOMMENT

