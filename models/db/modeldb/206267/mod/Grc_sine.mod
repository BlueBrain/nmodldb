TITLE Electrode for sinusoidal current clamping

COMMENT
	Author: A. Fontana
	Last revised: 28.3.99
ENDCOMMENT

NEURON {
	POINT_PROCESS sine
	RANGE del, dur, amp, i, freq, phasernd,amprnd,z
	ELECTRODE_CURRENT i
}
UNITS {
	(nA) = (nanoamp)
}

PARAMETER {
	PI = 3.141592
	del (ms)
	dur = 1e5 (ms)	
	amp = 0.0015 (nA)
	freq = 50 (1/ms)
}
ASSIGNED { 
	i 	(nA) 
	amprnd	(nA)
	phasernd
	z
}

INITIAL {
	i = 0
	z = 10
	while(z*z>9){ z = normrand(0,1) }
	phasernd=2*PI*scop_random()
	amprnd=amp*z
}

PROCEDURE seed(x) { set_seed(x) }

BREAKPOINT {
	at_time(del)
	at_time(del+dur)

	if (t < del + dur && t > del) {

		:printf("%g\n",poisrand(freq))

		i = amprnd*sin(2*PI*freq*t/1000+phasernd)	:*sin(2*PI*freq1*t/1000)


	}else{
		i = 0
	}
}
