COMMENT
Two state kinetic scheme synapse described by rise time tau1,
and decay time constant tau2. The normalized peak condunductance is 1.
Decay time MUST be greater than rise time.

The solution of A->G->bath with rate constants 1/tau1 and 1/tau2 is
 A = a*exp(-t/tau1) and
 G = a*tau2/(tau2-tau1)*(-exp(-t/tau1) + exp(-t/tau2))
	where tau1 < tau2

If tau2-tau1 -> 0 then we have a alphasynapse.
and if tau1 -> 0 then we have just single exponential decay.

The factor is evaluated in the
initial block such that an event of weight 1 generates a
peak conductance of 1.

Because the solution is a sum of exponentials, the
coupled equations can be solved as a pair of independent equations
by the more efficient cnexp method.

ENDCOMMENT

NEURON {
	POINT_PROCESS ScalExp2Syn
	: activation
	: --------------
	RANGE tau1, tau2, Ee, i
	RANGE Enumsyns			: number of synapses representd by this logsyn
	RANGE id					: the id number of the logsyn
	RANGE Eintrvl, NEproc	
	RANGE Eg
	POINTER proc_num
	: HSP
	: -----
	RANGE Egmax0, Egmax
	RANGE HSP_type			: 0-Vsoma; 1-local membrane potential
	RANGE Vtrg, V
	RANGE Etau, Eenable
	RANGE continuous_update, event_window : SF updated continuously or only within a window after a synaptic event
	RANGE order		: ODR order 0=additive; 1=multiplicative
	POINTER Vsoma
	: voltage averaging
	: -----------------
	RANGE vavg, avgstrt
	RANGE mavg, mavgstrt, mavgintrvl

	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
	(nS) = (nanosiemens)
}

PARAMETER {
	tau1=.1 (ms) <1e-9,1e9>
	tau2 = 10 (ms) <1e-9,1e9>
	Ee=0	(mV)
	Enumsyns = 1
	id = 0
	Eintrvl = 500 (ms)		: average interval between excitatory synaptic events
	NEproc				: number of processes from which an event is drawn every time step
	: HSP
	: -------
	Egmax0 = 1 (nS)
	Vtrg = -60 (mV)		: // target membrane potential
	Etau=1e5 (ms)			: // HSP time constant
	Eenable=1			: // enable HSP
	order=0				: // order of differential equation (order=0 => additive; order=1 => multiplicative)

	continuous_update=1		: // flag: 1=continuously update SF; 0=update only after synaptic event
	event_window=5 (ms)		: // window for SF update after a synaptic event
	HSP_type = 1
	
	: voltage averaging
	: -----------------
	avgstrt=0 (ms)			: // time to start averaging
	mavgstrt=0 (ms)			: // time to start moving average first interval
	mavgintrvl=5000 (ms)		: // size of interval for which average voltage is calculated
}

ASSIGNED {
	v (mV)
	i (nA)
	Eg (uS)
	factor
	Eorder
:	total (uS)
	vavg		  	: accumulated voltage average
	stp1			: number of simulation steps
	stp2			: number of simulation steps
	mavg			: moving average value at the end of an interval
	vintrvl			: sum of voltage from beginning of current interval to t
	n			: interval counter

	in_window		: flag: 1=within event window, update SF; 0=don't update SF
	synevent_time		: time of latest synaptic event
	proc_num
	Vsoma
	V
}

STATE {
	EA (uS)
	EB (uS)
	Egmax (nS)			: excitatory peack conductance
}

INITIAL {
	LOCAL j,tp
:	total = 0
	if (tau1/tau2 > .9999) { tau1 = .9999*tau2 }
	EA = 0
	EB = 0
	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
	factor = -exp(-tp/tau1) + exp(-tp/tau2)
	factor = 1/factor

	Egmax = Egmax0
	Eorder = Egmax^order
	if (HSP_type == 0) {
		V = Vsoma
	} else if (HSP_type == 1) {
		V = v
	}
	in_window = continuous_update	: intially, either in_window=1, or in_window=0 until a synaptic event occurs
	synevent_time = 0
	stp1 = 0
	stp2 = 0
	n = 1
	vavg = 0
	vintrvl = 0
	mavg = v
}

BREAKPOINT {
	Eg = EB - EA
	i = Eg*(v - Ee)
	SOLVE update METHOD after_cvode
	SOLVE state METHOD cnexp
}

DERIVATIVE state {
	EA' = -EA/tau1
	EB' = -EB/tau2
:	Egmax' = Eenable * in_window * Eorder * (Vtrg - V) / Etau
}

NET_RECEIVE(w) {
	LOCAL draw, j
	j = 0
	while (j < Enumsyns) {
		draw = floor(scop_random() * NEproc)	: draw a process number
		if (proc_num == draw) {
			synevent_time=t
			state_discontinuity(EA, EA + Egmax/1000*factor)
			state_discontinuity(EB, EB + Egmax/1000*factor)
			net_event(t)		
		}
		j = j + 1
	}
}

PROCEDURE update() {
	LOCAL sf
	: update Egmax
	: ---------
	if (Egmax<0) { Egmax=0 }
	Eorder = Egmax^order
	if (HSP_type == 0) {
		V = Vsoma
	} else if (HSP_type == 1) {
		V = v
	}
	sf = Vtrg - V
	if (Eenable) {
://		Egmax = Egmax * normrand(1, sf^2 / Etau)
		Egmax = Egmax * (1 + exprand(sf) / Etau)
	}
	: update EPSP_time
	: ----------------
	if (!continuous_update) {
		if (t-synevent_time<=event_window) { in_window=1
		} else { in_window=0 }
	}
	: update accumulated voltage average
	: ----------------------------------
	if (t>=avgstrt) {
   		vavg=(vavg*stp1+v)/(stp1+1)
		stp1=stp1+1
	}
	: update interval voltage average
	: -------------------------------
	if (t>=mavgstrt) {
      		vintrvl=vintrvl+v
		stp2=stp2+1
   		if (t>mavgstrt+n*mavgintrvl-dt/2) { : if the interval has ended (additional dt/2 otherwise t remains smaller than the expression)
      			if (mavgintrvl>0) { mavg=vintrvl/stp2 }
      			n=n+1
      			vintrvl=0
			stp2=0
   		}
	}
}
