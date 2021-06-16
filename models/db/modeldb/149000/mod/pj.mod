: $Id: netstim.mod 2212 2008-09-08 14:32:26Z hines $
: comments at end

NEURON	{ 
  ARTIFICIAL_CELL pj
  RANGE interval, number, start, count, intervalVect
  RANGE noise
  THREADSAFE : only true if every instance has its own distinct Random
  POINTER donotuse
}

PARAMETER {
	interval	= 15 (ms) <1e-9,1e9>: time between spikes (msec)
	number	    = 1 <0,1e9>	: number of spikes (independent of noise)
	start		= 1 (ms)	: start of first spike
	noise		= 0 <0,1>	: amount of randomness (0.0 - 1.0)
	meani       = 15 (ms)
	count       = 0
    intervalVect[100]
}

ASSIGNED {
	event (ms)
	on
	ispike
	donotuse
}

PROCEDURE seed(x) {
	set_seed(x)
}

INITIAL {
    count       = 0
      
    FROM i=0 TO 99 {
        intervalVect[i] = normrand(meani,5)    
    }
    
	on = 0 : off
	ispike = 0
	if (noise < 0) {
		noise = 0
	}
	if (noise > 1) {
		noise = 1
	}
	if (start >= 0 && number > 0) {
		on = 1
		: randomize the first spike so on average it occurs at
		: start + noise*interval
		event = start + invl(interval) - interval*(1. - noise)
		: but not earlier than 0
		if (event < 0) {
			event = 0
		}
		net_send(event, 3)
	}
}	

PROCEDURE init_sequence(t(ms)) {
	if (number > 0) {
		on = 1
		event = 0
		ispike = 0
	}
}

FUNCTION invl(mean (ms)) (ms) {
	if (mean <= 0.) {
		mean = .01 (ms) : I would worry if it were 0.
	}
	if (noise == 0) {
		invl = mean
	}else{
		invl = (1. - noise)*mean + noise*mean*erand()
	}
}
VERBATIM
double nrn_random_pick(void* r);
void* nrn_random_arg(int argpos);
ENDVERBATIM

FUNCTION erand() {
VERBATIM
	if (_p_donotuse) {
		/*
		:Supports separate independent but reproducible streams for
		: each instance. However, the corresponding hoc Random
		: distribution MUST be set to Random.negexp(1)
		*/
		_lerand = nrn_random_pick(_p_donotuse);
	}else{
		/* only can be used in main thread */
		if (_nt != nrn_threads) {
hoc_execerror("multithread random in NetStim"," only via hoc Random");
		}
ENDVERBATIM
		: the old standby. Cannot use if reproducible parallel sim
		: independent of nhost or which host this instance is on
		: is desired, since each instance on this cpu draws from
		: the same stream
		erand = exprand(1)
VERBATIM
	}
ENDVERBATIM
}

PROCEDURE noiseFromRandom() {
VERBATIM
 {
	void** pv = (void**)(&_p_donotuse);
	if (ifarg(1)) {
		*pv = nrn_random_arg(1);
	}else{
		*pv = (void*)0;
	}
 }
ENDVERBATIM
}

PROCEDURE next_invl() {
	if (number > 0) {
		event = invl(interval)
	}
	if (ispike >= number) {
		ispike = 0
		number=normrand(15,5)
		if (number<=0) {number=-number+1}
		event = event + normrand(30,10)
		interval = intervalVect[count]
		count = count + 1
		:interval = normrand(meani,5)
	}
}

NET_RECEIVE (w) {
	if (flag == 0) { : external event
		if (w > 0 && on == 0) { : turn on spike sequence
			: but not if a netsend is on the queue
			init_sequence(t)
			: randomize the first spike so on average it occurs at
			: noise*interval (most likely interval is always 0)
			next_invl()
			event = event - interval*(1. - noise)
			net_send(event, 1)
		}else if (w < 0) { : turn off spiking definitively
			on = 0
		}
	}
	if (flag == 3) { : from INITIAL
		if (on == 1) { : but ignore if turned off by external event
			init_sequence(t)
			net_send(0, 1)
		}
	}
	if (flag == 1 && on == 1) {
		ispike = ispike + 1
		net_event(t)
		next_invl()
		if (on == 1) {
			net_send(event, 1)
		}
	}
}
