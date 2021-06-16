TITLE GCa.mod

COMMENT

GCa.mod
conductance read from a file, starting after vthreshold is exceeded

ENDCOMMENT

NEURON {
	SUFFIX GCa
	USEION ca READ eca WRITE ica
	RANGE vthreshold, gcabar, gt, initflag, del
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
}

PARAMETER {
	vthreshold  (mV)
	gcabar (S/cm2)
}

ASSIGNED {
	v (mV)
	eca (mV)
	ica (mA/cm2)
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
		gt = stepforce(t - del, "gCa.dat")
		ica = gcabar*gt*(v - eca)
	}else{
		gt = 0
		ica = 0
	}
}

PROCEDURE check() {
	if (initflag == 0 && v > vthreshold) {
		initflag = 1
		del = t
	}
}
