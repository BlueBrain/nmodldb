COMMENT
Simple model of synaptic modulation of a second messenger "msg" by dopamine.
msg approaches its steady state msginf with a single time constant tau.
Initially msg has value msginf = 1, and tau is tau0
(default is 100 ms).

If an input event with positive weight w arrives at time t
  If msginf == 1
    msginf = 1+w
    tau = tau1
    send a self-event that will return at t+dur
  else
    move the self-event to t+dur (i.e. just prolong "on" phase)

If a self-event arrives
   msginf = 1
   tau = tau0

An input event with negative weight has no effect.


Affects other mechanisms via POINTER variables linked to msg.


This model has simple dynamics and no stream-specificity.
ENDCOMMENT

NEURON {
	POINT_PROCESS DAsyn
	RANGE tau, tau0, tau1, dur, msg, msginf
}

PARAMETER {
	tau0 = 100 (ms)
	tau1 = 70 (ms)
	dur = 600 (ms)
}

ASSIGNED {
	tau (ms)
	msginf (1)
}

STATE {
	msg (1)
}

INITIAL {
	tau = tau0
	msginf = 1
	msg = 1
}

BREAKPOINT {
	SOLVE state METHOD cnexp
}

DERIVATIVE state {
	msg' = (msginf - msg)/tau
}

NET_RECEIVE (w(1)) {
  if (flag == 0) { : this is an input event
	if (w>0) { : ignore events with nonpositive weight
		if (msginf==1) {
			msginf = 1 + w
			tau = tau1
			net_send(dur, 1)
		} else {
			net_move(t + dur)			
		}
	}
  } else { : this is a self-event
	msginf = 1
	tau = tau0
  }
}
