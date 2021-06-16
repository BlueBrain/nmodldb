
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS VRamp
	ELECTRODE_CURRENT i
	RANGE rs, vc, i
    RANGE del,vo1,vf1,dur1,vo2,vf2,dur2,vo3,vf3,dur3,vo4,vf4,dur4,vo5,vf5,dur5
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)

}


PARAMETER {
    del = 0 (ms)
    vo1 = 0 (mV)
    vf1 = 0 (mV)
    dur1 = 0 (ms)
    vo2 = 0 (mV)
    vf2 = 0 (mV)
    dur2 = 0 (ms)
    vo3 = 0 (mV)
    vf3 = 0 (mV)
    dur3 = 0 (ms)
    vo4 = 0 (mV)
    vf4 = 0 (mV)
    dur4 = 0 (ms)
    vo5 = 0 (mV)
    vf5 = 0 (mV)
    dur5 = 0 (ms)
    rs = 1e-07 (megohm) <1e-9, 1e9>
}

ASSIGNED {
	v (mV)	: automatically v + vext when extracellular is present
	i (nA)
	vc (mV)
	on
    t1 (ms)
    t2 (ms)
    t3 (ms)
    t4 (ms)
    t5 (ms) 
    t6 (ms)
    m1 (mV/ms)
    m2 (mV/ms)
    m3 (mV/ms)
    m4 (mV/ms)
    m5 (mV/ms)
}

STATE { }


INITIAL {
	on = 0
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

:PROCEDURE states() { }

PROCEDURE vstim() {
	if (dur1 > 0) { m1 = (vf1-vo1)/dur1 } else { m1 = 0 }
    if (dur2 > 0) { m2 = (vf2-vo2)/dur2 } else { m2 = 0 }
    if (dur3 > 0) { m3 = (vf3-vo3)/dur3 } else { m3 = 0 }
    if (dur4 > 0) { m4 = (vf4-vo4)/dur4 } else { m4 = 0 }
    if (dur5 > 0) { m5 = (vf5-vo5)/dur5 } else { m5 = 0 }
    on = 1
	t1 = del
    t2 = t1 + dur1
    t3 = t2 + dur2
    t4 = t3 + dur3
    t5 = t4 + dur4
    t6 = t5 + dur5
    at_time(t1)
    at_time(t2)
    at_time(t3)
    at_time(t4)
    at_time(t5)
    if (t>t1 && t<=t2) { vc = vo1 + m1*(t - t1) }
    if (t>t2 && t<=t3) { vc = vo2 + m2*(t - t2) }
    if (t>t3 && t<=t4) { vc = vo3 + m3*(t - t3) }
    if (t>t4 && t<=t5) { vc = vo4 + m4*(t - t4) }
	if (t>t5 && t<=t6) { vc = vo5 + m5*(t - t5) }
    icur()
}