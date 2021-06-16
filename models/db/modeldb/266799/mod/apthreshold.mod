NEURON {
	POINT_PROCESS APThreshold
	RANGE n, slopeth, ampth, time, vth, firing1, firing2
}

UNITS {
	(mV) = (millivolt)
}

PARAMETER {
	n
	slopeth = 10 (mV/ms)
	ampth = -20  (mV)
	time (ms)
	vth (mV)
	v (mV)
	dt (ms)
}

ASSIGNED {
	firing1
	firing2
	spacet
	spaceth
	v1 (mV)
}

VERBATIM
extern void vector_resize();
extern double* vector_vec();
extern void* vector_arg();
ENDVERBATIM

INITIAL {
	n = 0
	firing1 = 0
	firing2 = 0
	v1 = 0 (mV)
VERBATIM
	{ void* vvt; void* vvth;
		vvt = *((void**)(&spacet));
		vvth = *((void**)(&spaceth));
		if (vvt) {
			vector_resize(vvt, 0);
		}
		if (vvth) {
			vector_resize(vvth, 0);
		}
	}
ENDVERBATIM
	check()
}

BREAKPOINT {
	SOLVE check METHOD after_cvode
}

PROCEDURE check() {
VERBATIM
	if ((v - v1)/dt >= slopeth && !firing1) {
		firing1 = 1;
		time = t;
		vth = v;
	}
	if (v >= ampth && firing1 && !firing2) {
		int size; double* px; void* vvt; void* vvth;
		firing2 = 1;
		n += 1.;
		vvt = *((void**)(&spacet));
		vvth = *((void**)(&spaceth));
		if (vvt) {
			size = (int)n;
			vector_resize(vvt, size);
			px = vector_vec(vvt);
			px[size-1] = time;
		}
		if (vvth) {
			size = (int)n;
			vector_resize(vvth, size);
			px = vector_vec(vvth);
			px[size-1] = vth;
		}
	}
	if (firing1 && (v - v1)/dt < slopeth && t > time) {
		firing1 = 0;
	}
	if (firing2 && v < ampth && t > time) {
		firing2 = 0;
		firing1 = 0;
	}
	v1 = v;
ENDVERBATIM
}

PROCEDURE spiketimes() {
VERBATIM
	extern void* vector_arg();
	void** vvt;
	vvt = (void**)(&spacet);
	*vvt = (void*)0;
	if (ifarg(1)) {
		*vvt = vector_arg(1);
	}
ENDVERBATIM
}

PROCEDURE thresholds() {
VERBATIM
	extern void* vector_arg();
	void** vvth;
	vvth = (void**)(&spaceth);
	*vvth = (void*)0;
	if (ifarg(1)) {
		*vvth = vector_arg(1);
	}
ENDVERBATIM
}
