COMMENT
*************** M.Migliore 2001 ************
Pulses of in (0.2ms) are generated using a modified Poisson
distribution of arrival times with mean frequency drawn from a normal distr.
At this time it cannot be started from the INITIAL block using the variable time step. 
An external event must be used.
ENDCOMMENT
					       
NEURON {
	POINT_PROCESS pnoise
	RANGE e,onset,gavr,count,alp,eps,bet,rfreq,meant,type,sd,start,seed
	NONSPECIFIC_CURRENT i
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
}

PARAMETER {
	e=0	(mV)
	v	(mV)
	count=0
	onset=0 (ms)
	alp=0.8  (ms)
	bet=3  (ms)
	eps=1  (ms)
	on=0.02 (nA)
        gavr=0 (umho)
	meant=100 (ms) 
	rfreq=0 (ms)
	type=-1
	sd=0
	start=0
	seed=1
}

ASSIGNED {
	i (nA)
	in 
}

STATE {	c (microsiemens) <1.e-4> n (microsiemens) <1.e-4>}

INITIAL {
	c=0
	n=0
	in=0
	getpar(sd)
	getpar(type)
}

BREAKPOINT {
	SOLVE states METHOD sparse
	i = c*(v - e)
}

KINETIC states {
	~ in <-> n (1/alp, 0)
	~ n  <-> c (1/bet, 0)
	~ c   ->   (1/eps)
	}

PROCEDURE getpar(r) {
	if (r < 0) {
		if (r==-1) {onset=-rfreq*log(1-scop_random())
		} else {
		onset=rfreq
		}
		} else {
		while (rfreq<=0) {rfreq=normrand(meant,r)}
		}
}
NET_RECEIVE (w) {
	if (flag == 0) { : from external
			net_send(onset, 1)	:starts the process
			getpar(type)		:sets next onset
	}
	if (flag == 1) {
		in=gavr*on
		net_send(0, 0)
		net_send(.2, 2)
	}
	if (flag == 2) {
		in=0
		count=count+1
	}
}



