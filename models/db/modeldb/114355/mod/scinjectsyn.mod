COMMENT
Constant conductance injecting synapse
ENDCOMMENT

NEURON {
	POINT_PROCESS ScalInjectSyn
	: activation
	: --------------
	RANGE Ee, Eg, i
	: HSP
	: -----
	RANGE Egmax0, Egmax
	RANGE HSP_type			: 0-Vsoma; 1-local membrane potential
	RANGE Vtrg, V
	RANGE Etau, Eenable
	RANGE order		: ODR order 0=additive; 1=multiplicative
	POINTER Vsoma
	: voltage averaging
	: -----------------
	RANGE vavg, avgstrt
	RANGE mavg, mavgstrt, mavgintrvl
	: not used here, but referred to from hoc
	: -------------------------------------------------------
	RANGE tau1, tau2, Enumsyns, id, Eintrvl, NEproc	
	POINTER proc_num
	RANGE continuous_update, event_window

	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
	(nS) = (nanosiemens)
}

PARAMETER {
	Ee=0	(mV)
	: HSP
	: -------
	Egmax0 = 1 (nS)
	Vtrg = -60 (mV)		: target membrane potential
	Etau=1e5 (ms)			: HSP time constant
	Eenable=1			: enable HSP
	order=0				: order of differential equation (order=0 => additive; order=1 => multiplicative)
	HSP_type = 1
	: voltage averaging
	: -----------------
	avgstrt=0 (ms)			: time to start averaging
	mavgstrt=0 (ms)			: time to start moving average first interval
	mavgintrvl=5000 (ms)		: size of interval for which average voltage is calculated
	: not used but referred from hoc
	: -------------------------------------------
	tau1=.1 (ms) <1e-9,1e9>
	tau2 = 10 (ms) <1e-9,1e9>
	Enumsyns = 1
	id = 0
	Eintrvl = 500 (ms)
	NEproc
	continuous_update=1
	event_window=5 (ms)
}

ASSIGNED {
	v (mV)
	i (nA)
	Eg (uS)
	Eorder
	vavg		  	: accumulated voltage average
	stp1			: number of simulation steps
	stp2			: number of simulation steps
	mavg			: moving average value at the end of an interval
	vintrvl			: sum of voltage from beginning of current interval to t
	n			: interval counter
	Vsoma
	V
	: not used but referred from hoc
	: -------------------------------------------
	proc_num
}

STATE {
	Egmax (nS)			: excitatory peack conductance
}


INITIAL {
	Egmax = Egmax0
	Eorder = Egmax^order
	if (HSP_type == 0) {
		V = Vsoma
	} else if (HSP_type == 1) {
		V = v
	}
	stp1 = 0
	stp2 = 0
	n = 1
	vavg = 0
	vintrvl = 0
	mavg = v
}

BREAKPOINT {
	Eg = Egmax
	i = Eg*(v - Ee)
	SOLVE update METHOD after_cvode
	SOLVE state METHOD cnexp
}

DERIVATIVE state {
	Egmax' = Eenable * Eorder * (Vtrg - V) / Etau
}

NET_RECEIVE(w) {
	: not used, just needed for hoc recordings aritificial_vivo.hoc
}

PROCEDURE update() {
	: update Egmax
	: ---------
	if (Egmax<0) { Egmax=0 }
	Eorder = Egmax^order
	if (HSP_type == 0) {
		V = Vsoma
	} else if (HSP_type == 1) {
		V = v
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
