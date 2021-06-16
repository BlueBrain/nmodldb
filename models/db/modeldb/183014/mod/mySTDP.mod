: $Id: expsynstdp.mod,v 1.3 2012/02/17 15:45:45 samn Exp $ 
:
: basic STDP exponential synapse
: from http://www.neuron.yale.edu/neuron/static/news/stdp.mod
:

NEURON {
	POINT_PROCESS mySTDP : at a point
	RANGE tau, e, i, d, p, dtau, ptau : those values beong to this process and are not global
	NONSPECIFIC_CURRENT i : transmembrane current not necessarily calcium or specific ion
        GLOBAL verbose : all instances of this class can print you things if you put it to 1
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	tau = 0.1 (ms) <1e-9,1e9>
	e = 0	(mV)
	d = 0.2 <0,1>: depression factor (multiplicative to prevent < 0)
	p = 0.2 : potentiation factor (multiplicative)
	dtau = 34 (ms) : depression effectiveness time constant
	ptau = 17 (ms) : Bi & Poo (1998, 2001)
        verbose = 0 
}

ASSIGNED {
	v (mV)
	i (nA)
	tpost (ms)
}

STATE {
	g (uS)
}

INITIAL {
	g=0
	tpost = -1e9
	net_send(0, 1)
}

BREAKPOINT { : the value is updated at every time step, neuron automatically does the update for v so don't do it
	SOLVE state METHOD cnexp
	i = g*(v - e)
}

DERIVATIVE state { : put all your derivatives, they will be called at every time step
	g' = -g/tau
}

NET_RECEIVE(w (uS), A, tpre (ms)) { : presynaptic spike events are received and handled here, they always have flag 0
	INITIAL { A = 1  tpre = -1e9 }
	if (flag == 0) { : presynaptic spike  (after last post so depress)
          if(verbose) {printf("entry flag=%g t=%g w=%g A=%g tpre=%g tpost=%g\n", flag, t, w, A, tpre, tpost)}
		g = g + w*A
		tpre = t
		A = A * (1 - d*exp((tpost - t)/dtau))
	}else if (flag == 2) { : postsynaptic spike
          if(verbose) {printf("entry flag=%g t=%g tpost=%g\n", flag, t, tpost)}
		tpost = t : tracs the time of the post spike
		FOR_NETCONS(w1, A1, tp) { : also can hide NET_RECEIVE args, goes through all the sources, you can connect many spike sources to a single synapse
                  if(verbose) {printf("entry FOR_NETCONS w1=%g A1=%g tp=%g\n", w1, A1, tp)}
			A1 = A1*(1 + p*exp((tp - t)/ptau))
		}
	} else { : flag == 1 from INITIAL block
          if(verbose) {printf("entry flag=%g t=%g\n", flag, t)}
		WATCH (v > -20) 2 : if there is a postsynaptic spike, send an event with flag 2
	}
}
