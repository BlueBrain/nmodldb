: Na-channels from Aradi and Holmes 1999 as available in ModelDB (transformed into mod file)
: Note: "The Na and KDR parameters in the paper were chosen relative to a
: resting potential of 0 mV; we have adjusted them to correspond to a
: resting potential of -70 mV." 

NEURON {
	SUFFIX ichan3
	USEION na READ ena WRITE ina
	USEION k READ ek WRITE ik
	RANGE ik, ina, gnabar, gna,  gksbar, gks, gkabar, gka, gka,gks, gkfbar, gkf
	GLOBAL vshiftak,vshiftks,vshiftna
}

UNITS {
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER {
	gnabar (S/cm2) 
	gkfbar (S/cm2) 
	gksbar (S/cm2) 
	gkabar (S/cm2) 
	vshiftak = 0 (mV)
	vshiftks = 0 (mV)
	vshiftna = 0 (mV)
}

ASSIGNED {
	v (mV)
	ena (mV)
	ek (mV)
	ina (mA/cm2)
	ik (mA/cm2)
	gna (S/cm2)
	gkf (S/cm2)
	gks (S/cm2)
	gka (S/cm2)
}

STATE { m n1 h n2 k l }

BREAKPOINT {
SOLVE states METHOD cnexp
gna = gnabar * m^3*h
gkf = gkfbar * n1^4
gks = gksbar * n2^4
gka = gkabar * k * l
ina = gna * (v - ena)
ik = (gks+gka) * (v - ek) :gkf
}

INITIAL {
: Assume v has been constant for a long time
	m = am(v)/(am(v) + bm(v))
	h  = ah(v)/(ah(v) + bh(v))
	n1 = an1(v)/(an1(v) + bn1(v))
	n2 = an2(v)/(an2(v) + bn2(v))
	k = ak(v)/(ak(v) + bk(v))
	l = al(v)/(al(v) + bl(v))
}

DERIVATIVE states {
	: Computes state variables at present v & t
	m' = (1-m)*am(v) - m*bm(v)
	h' = (1-h)*ah(v) - h*bh(v)
	n1' = (1-n1)*an1(v) - n1*bn1(v)
	n2' = (1-n2)*an2(v) - n2*bn2(v)
	k' = (1-k)*ak(v) - k*bk(v)
	l' = (1-l)*al(v) - l*bl(v)
}

FUNCTION am(Vm (mV)) (/ms) {
	LOCAL x
	UNITSOFF
	x = 0.2*(Vm - (-45+vshiftna)) 
	if(fabs(x) > 1e-6) {
		am = 1.5*x/(1-exp(-x))
	}else{
		am = 1.5/(1+0.5*x)
	}
	UNITSON
}

FUNCTION bm(Vm (mV)) (/ms) {
	LOCAL x
	UNITSOFF
	x = -0.2*(Vm - (-17+vshiftna))
	if(fabs(x) > 1e-6) {
		bm = 1.5*x/(1-exp(-x))
	}else{
		bm = 1.5/(1+0.5*x)
	}
	UNITSON
}

FUNCTION ah(Vm (mV)) (/ms) {
	UNITSOFF
	ah = 0.23*exp(-0.05*(Vm - (-67+vshiftna)))
	UNITSON
}

FUNCTION bh(Vm (mV)) (/ms) {
	LOCAL x
	UNITSOFF
	x= -(-0.1)*((-14.5+vshiftna)-Vm)	
	bh = 3.33/(1+exp(x))
	UNITSON
}

FUNCTION an1(Vm (mV)) (/ms) {
	LOCAL x
	UNITSOFF
	x = 0.16667*(Vm+23)
	if(fabs(x) > 1e-6) {
		an1 = 0.42*x/(1 - exp(-x))
	}else{
		an1 = 0.42/(1 +0.5*x)
	}
	UNITSON
}

FUNCTION bn1(Vm (mV)) (/ms) {
	UNITSOFF
	bn1 = 0.264*exp(-0.025*(Vm+48)) 
	UNITSON
}

FUNCTION an2(Vm (mV)) (/ms) {
	LOCAL x
	UNITSOFF
	x = 0.16667*(Vm-(-35+vshiftks))
	if(fabs(x) > 1e-6) {
		an2 = 0.168*x/(1 - exp(-x))
	}else{
		an2 = 0.42/(1 +0.5*x)
	}
	UNITSON
}

FUNCTION bn2(Vm (mV)) (/ms) {
	UNITSOFF
	bn2 = 0.1056*exp(-0.025*(Vm-(-60+vshiftks)))
	UNITSON
}

FUNCTION ak(Vm (mV)) (/ms) {
	LOCAL x
	UNITSOFF
	x = 0.066667*(Vm-(-25+vshiftak))
	if(fabs(x) > 1e-6) {
		ak = 0.75*x/(1-exp(-x))
	}else{
		ak = 0.75/(1+0.5*x)
	}
	UNITSON
}

FUNCTION bk(Vm (mV)) (/ms) {
	LOCAL x
	UNITSOFF
	x = -0.125*(Vm-(-15 +vshiftak)) 
	if(fabs(x) > 1e-6) {
		bk = 0.8*x/(1-exp(-x))
	}else{
		bk = 0.8/(1+0.5*x)
	}
	UNITSON
}

FUNCTION al(Vm (mV)) (/ms) {
	UNITSOFF
	al = 0.00015*exp(-0.066667*(Vm -(-13+vshiftak)))
	UNITSON
}

FUNCTION bl(Vm (mV)) (/ms) {
	UNITSOFF
	bl = 0.06/(1+exp(-(-0.083333)*((-68+vshiftak)-Vm)))
	UNITSON
}