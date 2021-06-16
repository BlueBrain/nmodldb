TITLE vsource.mod
COMMENT
Patterned after svclmp.mod, a single electrode Voltage clamp mechanism.
Unlike svclmp.mod, this has one level--amp--and it delivers current i 
for t<=toff (even t<0).

Clamp is on at time 0, and off after time toff.
When clamp is off the injected current is 0.
i is the injected current, vc measures the control voltage.

For other comments and important implementational details, 
especially "why SOLVE icur METHOD after_cvode ???", 
see the comments in svclmp.mod
ENDCOMMENT

NEURON {
	POINT_PROCESS Vsource
	ELECTRODE_CURRENT i
	RANGE toff, amp, rs, vc, i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}


PARAMETER {
	rs = 1 (megohm) <1e-9, 1e9>
	toff (ms) 	  amp (mV)
}

ASSIGNED {
	v (mV)	: automatically v + vext when extracellular is present
	i (nA)
	vc (mV)
	on
}

INITIAL {
	on = 1
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

PROCEDURE vstim() {
	on = 1
	if (toff) {at_time(toff)}
	if (t < toff) {
		vc = amp
	}else {
		vc = 0
		on = 0
	}
	icur()
}
