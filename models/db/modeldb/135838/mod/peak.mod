TITLE peak.mod

COMMENT

pk: record peak time and peak value of membrane potential
Arnd Roth 25.9.1997
revised   10.9.2008

ENDCOMMENT

UNITS {
	(mV) = (millivolt)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX pk
	RANGE tmin, vmin
	RANGE tpeak, vpeak
	RANGE dvdt, d2vdt2, tdvdtpeak, dvdtpeak, td2vdt2peak, d2vdt2peak
}

PARAMETER {
	tmin	(ms)
	vmin	(mV)
	v	(mV)
}

ASSIGNED {
	tpeak	(ms)
	vpeak	(mV)
	v1	(mV)
	v2	(mV)
	v3	(mV)
	t1	(ms)
	t2	(ms)
	t3	(ms)
	initflag	(1)
	dvdt	(mV/ms)
	d2vdt2	(mV/ms2)
	tdvdtpeak	(ms)
	dvdtpeak	(mV/ms)
	td2vdt2peak	(ms)
	d2vdt2peak	(mV/ms2)
}

INITIAL {
	tpeak = 0	(ms)
	vpeak = -100	(mV)
	v1 = 0	(mV)
	v2 = 0	(mV)
	v3 = 0	(mV)
	t1 = 0	(mV)
	t2 = 0	(mV)
	t3 = 0	(mV)
	initflag = 0		(1)
	dvdt = 0		(mV/ms)
	d2vdt2 = 0		(mV/ms2)
	tdvdtpeak = 0		(ms)
	dvdtpeak = 0		(mV/ms)
	td2vdt2peak = 0		(ms)
	d2vdt2peak = 0		(mV/ms2)
	check()
}

BREAKPOINT {
	SOLVE check METHOD after_cvode
}

PROCEDURE check() {
	if (v > vpeak) {
		tpeak = t
		vpeak = v
	}
	v1 = v2
	v2 = v3
	v3 = v
	t1 = t2
	t2 = t3
	t3 = t
	if (t3 - t2 > 0) {
		dvdt = (v3 - v2)/(t3 - t2)
	}
	if (t3 - t2 > 0 && t2 - t1 > 0) {
		d2vdt2 = (2*(t3*(v2 - v1) + t2*(v1 - v3) + t1*(v3 - v2)))/((t1 - t2)*(t1 - t3)*(t2 - t3))
	}
	if (v > vmin && t > tmin) {
		if (dvdt > dvdtpeak) {
			tdvdtpeak = t
			dvdtpeak = dvdt
		}
		if (dvdt > 0 && d2vdt2 > d2vdt2peak) {
			td2vdt2peak = t
			d2vdt2peak = d2vdt2
		}
	}
}
