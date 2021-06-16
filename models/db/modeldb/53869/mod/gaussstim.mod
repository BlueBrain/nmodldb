COMMENT
This mod file is based on netstim.mod in NEURON
modified by Yi Zhou to generate Gaussian interval distributions
ENDCOMMENT


NEURON	{ 
  POINT_PROCESS GaussStim
  RANGE x
  RANGE interval, number, start,factor,refrac
  
  RANGE N_backward,N_forward,N_normal,N_total
  RANGE rnd
  RANGE rand
  
}

PARAMETER {
	interval	= 10 (ms) <1e-9,1e9>: time between spikes (msec)
	number	= 10000 <0,1e9>	: number of spikes
	start		= 10 (ms)	: start of first spike
	factor =4  : portion of std to mean
	refrac		=0.5 <0,3>    : absolute refractory period, up limit of freq=1kHz
							: refrac >=dt, otherwise x can't be reset to 0 
}

ASSIGNED {
	x
	event (ms)	
	on
	end (ms)
	m (ms)            : mean of Gaussian
	diff (ms)
	N_forward  : swap spike whose value exceed forwardly one interval 
	N_backward  : swap spike whose value exceed backwardly one interval 
	N_normal
	N_total
	rnd
	rand
}

PROCEDURE seed(x) {
	set_seed(x)
}

INITIAL {
	on = 0
	x = 0
	diff=0
	m=interval/2  : each T has normal distribution N(interval/2, interval/2/factor)
	N_forward=0
	N_backward=0
	N_normal=0
	N_total=0
	if (start >= 0 && number > 0) {
		: randomize the first spike so on average it occurs at start
		
		event = start + invl(m) 
		
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
		event = t
		end = t + 1e-6 + interval*(number-1)
		
	}
}

FUNCTION invl(mean (ms)) (ms) { LOCAL std
	if (mean <= 0.) {
		mean = .01 (ms) 
	}
	std=mean/factor  
	invl = normrand(mean, std)  : relative to current interval 
	
	if(invl>=interval) { 
		invl=fmod(invl,interval)
		N_forward=N_forward+1
		}else if(invl<0) { 

			invl=fmod(invl,interval)+interval
			
			N_backward=N_backward+1
			}else {
			N_normal=N_normal+1
			}
		
		diff=interval-invl
	
}

PROCEDURE event_time() {LOCAL diff2,T
        diff2=diff
	if (number > 0) {
	   T=invl(m)
	   rnd=T
	   if(T==0 && diff2==0) { T=T+dt } : previous and current spikes overlapped
	   
	   event = T+event + diff2    :compute absolute event time, relative to 0ms
	   
 	   N_total=N_total+1
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
		if(x == 0){  : after refractory
				rand=scop_random() : [0, 1] uniform
				x = rand
				net_event(t)
				event_time()  : after each spike call next spike time

		  		 if (event-t <= refrac+dt && event-t >= refrac) {
		      		 		event=event+dt    : happen when next event at reset edge time
		  		 }

		  		 if (on==1) {
					net_send(event - t, 1) : self-feed next event
		  				 }
				net_send(refrac, 2)  : refractory period
		
		} else if (x!=0) {
			net_event(t)
			event_time() : although this spike will be ignored, still call next spike : independent cycles 
			if (on == 1) {
				net_send(event - t, 1)
					  }
		}
		
	}
	if (flag == 2) {
		
		x = 0
		
	}
}


