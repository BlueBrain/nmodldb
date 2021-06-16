NEURON {
	POINT_PROCESS Randomizer
	RANGE Eintrvl			: mean excitatory ISI
	RANGE NEproc			: number of processes to generate
	POINTER proc_num		: current process number (global variable)
	RANGE special
}

PARAMETER {
	Eintrvl
	NEproc
	special = 0	: very non-elegant way to create the same input stream for 1 Hz
}

ASSIGNED {
	proc_num
	last_event
	tshift
	skip
	count
}

INITIAL {
	LOCAL j
	skip = 0
	count = 0
	if (special) {
		Eintrvl = Eintrvl / 2.5
	}
	last_event = 0
	tshift = 0
	j = 0
	while (j < NEproc) {
		net_send(exprand(Eintrvl),j)
		j = j + 1
	}
}

NET_RECEIVE(w) {
	if (t == last_event) {
		tshift = tshift + dt/10
		net_send(tshift,flag)
:		printf("Randomizer, event overlap: time = %g process %g\n",t,flag)
	} else {
		proc_num = flag
		net_send(exprand(Eintrvl),flag)
		: special case where subgroup has 1 Hz instead of 2 Hz and I want to create
		: the same input raster
		if (special) {
			count = count + 1
			if (count == 2 || count == 4) {
				skip = 0
			} else {
				skip = 1
			}
			if (count == 5) {
				count = 0
			}
		}
		if (!skip) {
			net_event(t)
		}
		tshift = 0
	}
	last_event = t
}

PROCEDURE seed(x) {
	set_seed(x)
}

