COMMENT
toggle a flag at beginning and end of a spike
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS Spike
	RANGE thresh
	RANGE spike
	RANGE spike_on
	RANGE spike_time
	RANGE spike_count
	RANGE spike_freq_isi	: a rough estimate based on 1/(last ISI)
	RANGE spike_freq_count	: a rough estimate based on spike_count/t
}

UNITS {
	(mV) = (millivolt)
}

PARAMETER {
	thresh = 0 (mV)	: threshold for spike identification
}

ASSIGNED {
	spike
	spike_on
	spike_time
	spike_freq_isi
	spike_freq_count
	spike_count
}

INITIAL {
	spike=0
	spike_time=0
	spike_freq_isi=0
	spike_on = 0
	spike_count=0
	spike_freq_count=0
	if (v > thresh) {
		spike_on = 1
	}
}

BREAKPOINT {
	SOLVE check
}

PROCEDURE check() {
	if (spike) {
		spike=0 : spike=1 just for the initial dt of the spike
	}
	if (spike_on && v < thresh) {
		spike_on = 0
	}
	if (!spike_on && v > thresh) {
		spike=1
		spike_on = 1
		spike_count=spike_count+1
		if (t) {
			spike_freq_isi=1000/(t-spike_time)
			spike_freq_count=spike_count/t*1000
		}
		spike_time = t
	}
}
