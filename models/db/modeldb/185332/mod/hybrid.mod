: Hybrid clamp point mechanism

TITLE hybrid.mod
COMMENT

ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS hybrid
	ELECTRODE_CURRENT i
	RANGE delay, dur, vc_amp, tot_dur, rs, vc, i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}


PARAMETER {
	rs = 1 (megohm) <1e-9, 1e9>
    delay = 4000 (ms)
    dur = 4000 (ms)
    vc_amp = -70 (mV)
    tot_dur = 25000 (ms)

}

ASSIGNED {
	v (mV)	: automatically v + vext when extracellular is present
	i (nA)
	vc (mV)
	tc2 (ms)
	tc3 (ms)
	on
}

INITIAL {
	
	on = 0
    i = 0
}

BREAKPOINT {
	SOLVE icur METHOD after_cvode
	vstim()
}

PROCEDURE icur() {
	if (on) {
		i = (vc - v)/rs
	}else{
		i = 0
	}
}

COMMENT
The SOLVE of icur() in the BREAKPOINT block is necessary to compute
i=(vc - v(t))/rs instead of i=(vc - v(t-dt))/rs
This is important for time varying vc because the actual i used in
the implicit method is equivalent to (vc - v(t)/rs due to the
calculation of di/dv from the BREAKPOINT block.
The reason this works is because the SOLVE statement in the BREAKPOINT block
is executed after the membrane potential is advanced.

It is a shame that vstim has to be called twice but putting the call
in a SOLVE block would cause playing a Vector into vc to be off by one
time step.
ENDCOMMENT

PROCEDURE vstim() {
    at_time(delay)
    at_time(delay+dur)
    
	if (t < delay) {
		on = 1
        vc = vc_amp
	}else if (t < delay+dur) {
        on = 0
        vc = 0
	}else {
		vc = vc_amp
		on = 1
	}
	icur()
}

