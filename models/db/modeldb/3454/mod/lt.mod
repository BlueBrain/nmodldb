: 	Low threshold potassium chanel from 
:	Contribution of the Kv3.1 potassium channel to high-frequency firing in mouse auditory neurones
:	Lu-Yang Wang, Li Gan, Ian D. Forsythe and Leonard K. Kaczmarek
:	J. Physiol (1998), 501.9, pp. 183-194

NEURON {
	SUFFIX LT
	USEION k READ ek WRITE ik
	RANGE gbar, g, ik
	GLOBAL linf, ltau, rinf, rtau, al, bl, ar, br
}

: area in paper is 1000 (um2) so divide our density parameters by 10

UNITS {
	(mV) = (millivolt)
	(S) = (mho)
	(mA) = (milliamp)
}

PARAMETER {
	gbar = .002 (S/cm2) : .02 (uS)
	gamma = .1

	kal = 1.2 (/ms)
	eal = .03512 (/mV)
	kbl = .2248 (/ms)
	ebl = -.0319 (/mV)

	kar = .0438 (/ms)
	ear = -.0053 (/mV)
	kbr = .0562 (/ms)
	ebr = -.0047 (/mV)
}

ASSIGNED {
	v (mV)
	ek (mV)
	ik (mA/cm2)

	linf
	ltau (ms)
	rinf
	rtau (ms)

	al (/ms)
	bl (/ms)
	ar (/ms)
	br (/ms)
}

STATE {
	l r
}

INITIAL {
	rates(v)
	l = linf
	r = rinf
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	ik = gbar*l*r*(v - ek) : pemdas may be a problem
}

DERIVATIVE state {
	rates(v)
	l' = (linf - l)/ltau
	r' = (rinf - r)/rtau
}

PROCEDURE rates(v(mV)) {
	al = kal*exp(eal*v)
	bl = kbl*exp(ebl*v)

	ar = kar*exp(ear*v)
	br = kbr*exp(ebr*v)

	linf = al/(al + bl)
	ltau = 1/(al + bl)
	rinf = ar/(ar + br)
	rtau = 1/(ar + br)
}

