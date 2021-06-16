
: Point process to record Cai amplitudes

NEURON {
	POINT_PROCESS CalAmp
	USEION ca READ cai	
	RANGE n, thresh, time, firing, max, high
}

UNITS {
 	(molar) = (1/liter)
  	(uM)    = (micromolar)
  	(mM)    = (millimolar)
}

PARAMETER {
	n
	thresh = 100e-6 (mM)
	time (ms)
}

ASSIGNED {
	firing
	space
	high
	max
	cai (mM)
}

VERBATIM
extern void vector_resize();
extern double* vector_vec();
extern void* vector_arg();
ENDVERBATIM

INITIAL {
	n = 0
	firing = 0
	high=0
VERBATIM
	{ void* vv;
		vv = *((void**)(&space));
		if (vv) {
			vector_resize(vv, 0);
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
	int size; double* px; void* vv;
	if (cai >= thresh && !firing) {
		firing = 1;
		time = t;
		high = 1;
		max=cai;
	}

	if(high) {
		if (cai<=thresh && t>time){
			n += 1.;
			vv = *((void**)(&space));
			if (vv) {
				size = (int)n;
				vector_resize(vv, size);
				px = vector_vec(vv);
				px[size-1] = max;
			}
			high=0;
		}	

		if(cai>max){
			max=cai;
		}	
	}		

	if (firing && cai < thresh && t > time) {
		firing = 0;
	}
ENDVERBATIM
}

PROCEDURE record() {
VERBATIM
	extern void* vector_arg();
	void** vv;
	vv = (void**)(&space);
	*vv = (void*)0;
	if (ifarg(1)) {
		*vv = vector_arg(1);
	}
ENDVERBATIM
}
