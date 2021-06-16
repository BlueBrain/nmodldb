TITLE GNa.mod

COMMENT

GNa.mod
conductance read from a file, starting after vthreshold is exceeded

ENDCOMMENT

NEURON {
	SUFFIX GNa
	USEION na READ ena WRITE ina
	RANGE vthreshold, gnabar, gt, initflag, del
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
}

PARAMETER {
	vthreshold  (mV)
	gnabar (S/cm2)
}

ASSIGNED {
	v (mV)
	ena (mV)
	ina (mA/cm2)
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
		gt = stepforce(t - del, "gNa.dat")
		ina = gnabar*gt*(v - ena)
	}else{
		gt = 0
		ina = 0
	}
}

PROCEDURE check() {
	if (initflag == 0 && v > vthreshold) {
		initflag = 1
		del = t
	}
}
