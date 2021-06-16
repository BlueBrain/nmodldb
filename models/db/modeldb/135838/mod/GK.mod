TITLE GK.mod

COMMENT

GK.mod
conductance read from a file, starting after vthreshold is exceeded

ENDCOMMENT

NEURON {
	SUFFIX GK
	USEION k READ ek WRITE ik
	RANGE vthreshold, gkbar, gt, initflag, del
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
}

PARAMETER {
	vthreshold  (mV)
	gkbar (S/cm2)
}

ASSIGNED {
	v (mV)
	ek (mV)
	ik (mA/cm2)
	initflag (1)
	del (ms)
	gt (1)
}

INITIAL {
	initflag = 0
	del = 0
	gt = 0
}

BREAKPOINT {
	SOLVE check METHOD after_cvode
	at_time(del)

	if (initflag) {
		gt = stepforce(t - del, "gK.dat")
		ik = gkbar*gt*(v - ek)
	}else{
		gt = 0
		ik = 0
	}
}

PROCEDURE check() {
	if (initflag == 0 && v > vthreshold) {
		initflag = 1
		del = t
	}
}
